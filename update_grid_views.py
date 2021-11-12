import os

import mobie
import pandas as pd


def fix_grid_table(table_folder):
    tab_path = os.path.join(table_folder, "default.tsv")
    default = pd.read_csv(tab_path, sep="\t")
    morpho = pd.read_csv(os.path.join(table_folder, "morphology.tsv"), sep="\t")
    default = default.merge(morpho)
    default = default.rename(columns={"grid_id": "annotation_id"})
    default.to_csv(tab_path, index=False, sep="\t")


def update_view(ds_name):
    cell_ids = (1, 2, 3)

    have_pyr, have_starch = False, False
    if ds_name in ("micromonas", "phaeodactylum"):
        have_pyr = True
    if ds_name == "micromonas":
        have_starch = True

    def _get_sources(cell_id):
        sources = [f"cell{cell_id}-raw", f"cell{cell_id}-mask", f"cell{cell_id}-organelles"]
        if have_pyr:
            sources += [f"cell{cell_id}-pyrenoid-mask"]
        if have_starch:
            sources += [f"cell{cell_id}-starch-mask"]
        return sources

    sources = [_get_sources(cell_id) for cell_id in cell_ids]
    display_groups = {f"cell{cell_id}-raw": "raw" for cell_id in cell_ids}
    display_groups.update(
        {f"cell{cell_id}-mask": "cell-masks" for cell_id in cell_ids}
    )
    display_groups.update(
        {f"cell{cell_id}-organelles": "organelles" for cell_id in cell_ids}
    )
    if have_pyr:
        display_groups.update(
            {f"cell{cell_id}-pyrenoid-mask": "pyrenoid-masks" for cell_id in cell_ids}
        )
    if have_starch:
        display_groups.update(
            {f"cell{cell_id}-starch-mask": "starch-masks" for cell_id in cell_ids}
        )

    table_folder = os.path.join("data", ds_name, "tables", "all-cells")
    fix_grid_table(table_folder)

    mobie.metadata.add_grid_bookmark(
        dataset_folder=f"./data/{ds_name}",
        name="all-cells",
        sources=sources,
        display_groups=display_groups,
        overwrite=True
    )


def update_grid_view(ds_folder):
    ds_name = os.path.split(ds_folder)[1]
    metadata = mobie.metadata.read_dataset_metadata(ds_folder)
    views = metadata["views"]
    if "all-cells" in views:
        update_view(ds_name)


def update_all_grid_views():
    datasets = mobie.metadata.read_project_metadata("./data")["datasets"]
    for ds in datasets:
        print(ds)
        update_grid_view("./data/" + ds)


update_all_grid_views()
