import os
import numpy as np
import pandas as pd


def write_merged_dataset_tables(project_folder_path, table):
    for dataset in ["emiliania", "nanochloropsis", "pelagomonas", "symbiodinium"]:
        cut_table = table.loc[table["species"].str.lower() == dataset, :].drop(columns="species")
        cut_table = cut_table.rename(columns={"cell_id": "label_id"})
        cut_table["label_id"] = cut_table["label_id"].astype(float)

        cell_mask_dir_path = os.path.join(project_folder_path, "data", dataset, "tables", "cell-mask")
        morphology_table_path = os.path.join(cell_mask_dir_path, "morphology.tsv")
        cut_table.to_csv(morphology_table_path, sep="\t", index=False)


def write_unmerged_dataset_tables(project_folder_path, table):
    for dataset in ["galdieria", "micromonas", "phaeodactylum"]:
        cut_table = table.loc[table["species"].str.lower() == dataset, :].drop(columns="species")
        cut_table = cut_table.rename(columns={"cell_id": "cell"})

        grid_dir_path = os.path.join(project_folder_path, "data", dataset, "tables", "all-cells")
        default_table_path = os.path.join(grid_dir_path, "default.tsv")
        morphology_table_path = os.path.join(grid_dir_path, "morphology.tsv")

        default_table = pd.read_csv(default_table_path, sep="\t")
        merged_table = default_table.merge(cut_table, how="left", on="cell", validate="one_to_one")

        merged_table = merged_table[["grid_id", "volume_micron_cubed", "surface_area_micron_squared"]]
        merged_table.to_csv(morphology_table_path, sep="\t", index=False)


def write_morphology_tables(project_folder_path, table_path):
    table = pd.read_csv(table_path, sep="\t")
    write_merged_dataset_tables(project_folder_path, table)
    write_unmerged_dataset_tables(project_folder_path, table)


if __name__ == '__main__':
    write_morphology_tables('.', './misc/collated_plankton_cell_data.tsv')
