from copy import deepcopy
import mobie


def update_grid_view(ds_folder):
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    views = metadata['views']
    new_views = {}

    for name, view in views.items():
        if 'sourceTransforms' in view:
            new_view = deepcopy(view)
            new_view['sourceTransforms'][0]['grid']['tables'].append('morphology.tsv')
            new_views[name] = new_view
        else:
            new_views[name] = view

    metadata['views'] = new_views
    mobie.metadata.write_dataset_metadata(ds_folder, metadata)


def update_all_grid_views():
    datasets = mobie.metadata.read_project_metadata('./data')['datasets']
    for ds in datasets:
        print(ds)
        update_grid_view('./data/' + ds)


update_all_grid_views()
