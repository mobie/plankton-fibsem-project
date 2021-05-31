import os
import h5py
import numpy as np
import mobie
import nrrd

root = '/g/schwab/Kimberly/Publications/MoBIE_paper/Clarisse_project/raw'
scale_factors = 5 * [[2, 2, 2]]
chunks = (64,) * 3


def copy_data(path, ds_name, name):
    data, header = nrrd.read(path, index_order='C')
    vrange = (data.min(), data.max())
    print(vrange)
    resolution = np.diag(header['space directions'])[::-1] / 1000.
    tmp_path = f'./tmp_data/{ds_name}'
    os.makedirs(tmp_path, exist_ok=True)
    tmp_path = os.path.join(tmp_path, name + '.h5')
    if os.path.exists(tmp_path):
        return tmp_path, resolution.tolist(), vrange
    with h5py.File(tmp_path, 'w') as f:
        f.create_dataset('data', data=data, chunks=chunks)
    return tmp_path, resolution.tolist(), vrange


def add_data(path, ds_name, name, menu):
    tmp_path, resolution, vrange = copy_data(path, ds_name, name)
    mobie.add_image(tmp_path, "data", "./data", ds_name,
                    image_name=name, resolution=resolution,
                    scale_factors=scale_factors, chunks=chunks,
                    menu_name=menu, target='local',
                    max_jobs=16)


def create_dataset(name):
    ds_name = name.split('-')[0].lower()
    # for cell_id in (1, 2, 3):
    for cell_id in (1,):

        # add the em data
        raw_path = os.path.join(root, name, f'Data{cell_id}/Data{cell_id}.nrrd')
        assert os.path.exists(raw_path), raw_path
        raw_name = f'raw-cell{cell_id}'
        add_data(raw_path, ds_name, raw_name, 'fibsem-raw')

        # add the cell mask
        seg_path = os.path.join(root, name, f'Data{cell_id}/Segmented-files/Cell{cell_id}.nrrd')
        seg_name = f'cell-mask-cell{cell_id}'
        add_data(seg_path, ds_name, seg_name, 'fibsem-segmentation')

        # add the organelle segmentations
        organelle_paths = [
            f'Cell-coccolith{cell_id}.nrrd',
            f'Mitochondria{cell_id}.nrrd',
            f'Nucleus{cell_id}.nrrd',
            f'Plastid{cell_id}.nrrd',

        ]
        organelle_names = ['coccolith', 'mitochondria', 'nucleus', 'plastid']


def create_project():
    create_dataset('Emiliania-cell')


if __name__ == '__main__':
    create_project()
