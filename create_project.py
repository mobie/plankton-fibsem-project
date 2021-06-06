import os
from glob import glob

import mobie
import nrrd
import numpy as np
import pandas as pd
import z5py

root = '/g/schwab/Kimberly/Publications/MoBIE_paper/Clarisse_project/raw'
scale_factors = 4 * [[2, 2, 2]]
chunks = (64, 64, 64)


#
# utils
#


def copy_data(path, ds_name, name, invert=False, binarize=False):
    data, header = nrrd.read(path, index_order='C')
    if invert:
        data = data.max() - data
    if binarize:
        data = (data > 0).astype(data.dtype)
    vrange = [data.min(), data.max()]
    resolution = np.diag(header['space directions'])[::-1] / 1000.
    tmp_path = f'./tmp_data/{ds_name}'
    os.makedirs(tmp_path, exist_ok=True)
    tmp_path = os.path.join(tmp_path, name + '.n5')
    if os.path.exists(tmp_path):
        return tmp_path, resolution.tolist(), vrange
    with z5py.File(tmp_path, 'w') as f:
        f.create_dataset('data', data=data, chunks=chunks, n_threads=8)
    return tmp_path, resolution.tolist(), vrange


def add_data(path, ds_name, name, menu, invert=False, binarize=False, opacity=1.0):
    tmp_path, resolution, vrange = copy_data(path, ds_name, name, invert=invert, binarize=binarize)
    view = mobie.metadata.get_default_view('image', name, contrastLimits=vrange, opacity=opacity)
    mobie.add_image(tmp_path, "data", "./data", ds_name,
                    image_name=name, resolution=resolution,
                    scale_factors=scale_factors, chunks=chunks,
                    menu_name=menu, target='local',
                    max_jobs=16, view=view)


def get_organelles(name, cell_id):
    this_folder = os.path.join(root, name, f'Data{cell_id}/Segmented-files')
    all_paths = glob(os.path.join(this_folder, '*.nrrd'))

    organelle_paths = []
    organelle_names = []
    for path in all_paths:
        name = os.path.splitext(os.path.split(path)[1])[0]
        if name.startswith('Cell') or name.startswith('Starch') or name.startswith('Pyrenoid'):
            continue
        else:
            organelle_paths.append(path)
            organelle_names.append(name[:-1].lower())
    return organelle_paths, organelle_names


#
# split dataset:
# these use different raw data for the 3 cells with segmentation
#


def copy_segmentations(paths, ds_name, seg_name):
    data, header = nrrd.read(paths[0], index_order='C')
    data = (data > 0).astype('uint32')
    resolution = np.diag(header['space directions'])[::-1] / 1000.
    for ii, path in enumerate(paths[1:], 2):
        this_data, _ = nrrd.read(path, index_order='C')
        assert this_data.shape == data.shape
        data[this_data > 0] = ii
    assert np.array_equal(np.unique(data), np.arange(0, len(paths) + 1)),\
        f"{np.unique(data)}, {np.arange(0, len(paths) + 1)}"
    tmp_path = f'./tmp_data/{ds_name}'
    os.makedirs(tmp_path, exist_ok=True)
    tmp_path = os.path.join(tmp_path, seg_name + '.n5')
    if os.path.exists(tmp_path):
        return tmp_path, resolution.tolist()
    with z5py.File(tmp_path, 'w') as f:
        f.create_dataset('data', data=data, chunks=chunks, n_threads=8)
    return tmp_path, resolution.tolist()


def add_segmentations(paths, class_names, ds_name, seg_name, menu, opacity=0.5):
    assert len(paths) == len(class_names)
    tmp_path, resolution = copy_segmentations(paths,  ds_name, seg_name)
    view = mobie.metadata.get_default_view('segmentation', seg_name, opacity=opacity)
    mobie.add_segmentation(tmp_path, 'data', './data', ds_name,
                           segmentation_name=seg_name, resolution=resolution,
                           scale_factors=scale_factors, chunks=chunks,
                           menu_name=menu, target='local', max_jobs=16,
                           view=view)
    # add semantic column
    table_path = f'./data/{ds_name}/tables/{seg_name}/default.tsv'
    assert os.path.exists(table_path)
    table = pd.read_csv(table_path)
    table["organelle"] = class_names
    table.to_csv(table_path, sep='\t', index=False)


def _add_grid_view(ds_name, have_pyr, have_starch):
    cell_ids = (1, 2, 3)

    def _get_sources(cell_id):
        sources = [f'cell{cell_id}-raw', f'cell{cell_id}-mask', f'cell{cell_id}-organelles']
        if have_pyr:
            sources += [f'cell{cell_id}-pyrenoid-mask']
        if have_starch:
            sources += [f'cell{cell_id}-starch-mask']
        return sources

    sources = [_get_sources(cell_id) for cell_id in cell_ids]
    display_groups = {f'cell{cell_id}-raw': 'raw' for cell_id in cell_ids}
    display_groups.update(
        {f'cell{cell_id}-mask': 'cell-masks' for cell_id in cell_ids}
    )
    display_groups.update(
        {f'cell{cell_id}-organelles': 'organelles' for cell_id in cell_ids}
    )
    if have_pyr:
        display_groups.update(
            {f'cell{cell_id}-pyrenoid-mask': 'pyrenoid-masks' for cell_id in cell_ids}
        )
    if have_starch:
        display_groups.update(
            {f'cell{cell_id}-starch-mask': 'starch-masks' for cell_id in cell_ids}
        )
    mobie.metadata.add_grid_bookmark(
        dataset_folder=f'./data/{ds_name}',
        name='all-cells',
        sources=sources,
        display_groups=display_groups,
        overwrite=True
    )

    # fix the grid view table
    table_path = os.path.join('data', ds_name, 'tables', 'all-cells', 'default.tsv')
    table = pd.read_csv(table_path, sep='\t')
    table['cell'] = np.arange(1, len(table) + 1)
    table.to_csv(table_path, sep='\t', index=False)


def create_split_dataset(name):
    ds_name = name.split('-')[0].lower()
    if mobie.metadata.dataset_exists('./data', ds_name):
        print("Skipping", name)
        return

    have_pyr, have_starch = False, False
    for cell_id in (1, 2, 3):

        # add the em data
        raw_path = os.path.join(root, name, f'Data{cell_id}/Data{cell_id}.nrrd')
        assert os.path.exists(raw_path), raw_path
        raw_name = f'cell{cell_id}-raw'
        add_data(raw_path, ds_name, raw_name, menu='fibsem', invert=True)

        # add the cell mask
        seg_path = os.path.join(root, name, f'Data{cell_id}/Segmented-files/Cell{cell_id}.nrrd')
        seg_name = f'cell{cell_id}-mask'
        add_data(seg_path, ds_name, seg_name, menu='segmentation', binarize=True, opacity=0.5)

        # add additional masks depending on the dataset
        if name in ('Micromonas-cell', 'Phaeodactylum-cell'):
            seg_path = os.path.join(root, name, f'Data{cell_id}/Segmented-files/Pyrenoid{cell_id}.nrrd')
            seg_name = f'cell{cell_id}-pyrenoid-mask'
            add_data(seg_path, ds_name, seg_name, menu='segmentation', binarize=True, opacity=0.5)
            have_pyr = True
        if name == 'Micromonas-cell':
            seg_path = os.path.join(root, name, f'Data{cell_id}/Segmented-files/Starch{cell_id}.nrrd')
            seg_name = f'cell{cell_id}-starch-mask'
            add_data(seg_path, ds_name, seg_name, menu='segmentation', binarize=True, opacity=0.5)
            have_starch = True

        # add the organelle segmentations
        organelle_paths, organelle_names = get_organelles(name, cell_id)
        assert len(organelle_names) == len(organelle_paths)
        seg_name = f'cell{cell_id}-organelles'
        add_segmentations(organelle_paths, organelle_names, ds_name, seg_name, menu='segmentation', opacity=0.25)

    _add_grid_view(ds_name, have_pyr, have_starch)


def create_split_project():
    split_sets = [
        'Galdieria-cell',
        'Micromonas-cell',
        'Phaeodactylum-cell'
    ]
    for split_set in split_sets:
        create_split_dataset(split_set)


#
# merged dataset:
# these use the same raw data for the three cells with segmentations
#


def check_shapes(name):
    root_ds = os.path.join(root, name)
    shape = None
    for cell_id in (1, 2, 3):
        path = os.path.join(root_ds, f'Data{cell_id}/Data{cell_id}.nrrd')
        if not os.path.exists(path):
            path = os.path.join(root, name, f'Data{cell_id}/Cell{cell_id}.nrrd')
        assert os.path.exists(path)

        this_shape = nrrd.read_header(path)['sizes']
        if shape is None:
            shape = this_shape
        else:
            assert np.array_equal(shape, this_shape)


def add_merged_cell_mask(name, ds_name, semantic_name='Cell', seg_name='cell-mask'):
    root_ds = os.path.join(root, name)
    data = None
    for cell_id in (1, 2, 3):
        path = os.path.join(root_ds, f'Data{cell_id}/Segmented-files/{semantic_name}{cell_id}.nrrd')
        this_data, header = nrrd.read(path, index_order='C')
        resolution = (np.diag(header['space directions'])[::-1] / 1000.).tolist()

        if data is None:
            data = (this_data > 0).astype('uint32') * cell_id
        else:
            data[this_data > 0] = cell_id

    tmp_path = f'./tmp_data/{name}/merged_{seg_name}.n5'
    with z5py.File(tmp_path, 'w') as f:
        f.create_dataset('data', data=data, compression='gzip', chunks=chunks, n_threads=8)
    mobie.add_segmentation(tmp_path, 'data', './data', ds_name,
                           segmentation_name=seg_name, resolution=resolution,
                           scale_factors=scale_factors, chunks=chunks,
                           menu_name="segmentation", target='local', max_jobs=16)


def add_merged_organelles(name, ds_name):
    data = None
    organelle_names = None

    current_id = 1
    ids_to_cells, ids_to_organelles = [], []
    for cell_id in (1, 2, 3):
        this_paths, this_names = get_organelles(name, cell_id)
        if organelle_names is None:
            organelle_names = this_names
        else:
            assert all(n1 == n2 for n1, n2 in zip(organelle_names, this_names))

        for path, organelle_name in zip(this_paths, organelle_names):
            this_data, header = nrrd.read(path, index_order='C')
            resolution = (np.diag(header['space directions'])[::-1] / 1000.).tolist()
            if data is None:
                data = (this_data > 0).astype('uint32') * current_id
            else:
                data[this_data > 0] = current_id

            ids_to_cells.append(cell_id)
            ids_to_organelles.append(organelle_name)
            current_id += 1

    tmp_path = f'./tmp_data/{name}/merged_organelles.n5'
    with z5py.File(tmp_path, 'w') as f:
        f.create_dataset('data', data=data, compression='gzip', chunks=chunks, n_threads=8)
    seg_name = "organelles"
    mobie.add_segmentation(tmp_path, 'data', './data', ds_name,
                           segmentation_name=seg_name, resolution=resolution,
                           scale_factors=scale_factors, chunks=chunks,
                           menu_name="segmentation", target='local', max_jobs=16)

    # add columns
    table_path = f'./data/{ds_name}/tables/{seg_name}/default.tsv'
    assert os.path.exists(table_path)
    table = pd.read_csv(table_path)
    table["cell_id"] = ids_to_cells
    table["organelle"] = ids_to_organelles
    table.to_csv(table_path, sep='\t', index=False)


def create_merged_dataset(name):
    ds_name = name.split('-')[0].lower()
    if mobie.metadata.dataset_exists('./data', ds_name):
        print("Skipping", name)
        return
    check_shapes(name)

    raw_path = os.path.join(root, name, 'Data1/Data1.nrrd')
    if not os.path.exists(raw_path):
        raw_path = os.path.join(root, name, 'Data1/Cell1.nrrd')
    assert os.path.exists(raw_path)
    add_data(raw_path, ds_name, name='raw', menu='fibsem', invert=True)

    add_merged_cell_mask(name, ds_name)
    # add the cocolith mask (only for Emiliania)
    if name == 'Emiliania-cell':
        add_merged_cell_mask(name, ds_name, semantic_name='Cell-coccolith', seg_name='cocolith-mask')

    add_merged_organelles(name, ds_name)


def create_merged_project():
    merged_sets = [
        'Emiliania-cell',
        'Nanochloropsis-cell',
        'Pelagomonas-cell',
        'Symbiodinium-cell'
    ]
    for merged_set in merged_sets:
        create_merged_dataset(merged_set)


# TODO add metadata for the remote s3 store and upload it
if __name__ == '__main__':
    # create_merged_project()
    create_split_project()
