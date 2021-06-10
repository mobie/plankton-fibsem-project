import os
import numpy as np
import pandas as pd


def write_merged_dataset_tables(project_folder_path, table):
    for dataset in ["emiliania", "nanochloropsis", "pelagomonas", "symbiodinium"]:
        cut_table = table.loc[table["species"].str.lower() == dataset, :].drop(columns="species")

        organelles_dir_path = os.path.join(project_folder_path, "data", dataset, "tables", "organelles")
        default_table_path = os.path.join(organelles_dir_path, "default.tsv")
        morphology_table_path = os.path.join(organelles_dir_path, "morphology.tsv")

        default_table = pd.read_csv(default_table_path, sep="\t")
        merged_table = default_table.merge(cut_table, how="left", on=["organelle", "cell_id"], validate="one_to_one")

        merged_table = merged_table[["label_id", "volume_micron_cubed", "surface_area_micron_squared", "ratio_percent"]]
        merged_table.to_csv(morphology_table_path, sep="\t", index=False)


def write_unmerged_dataset_tables(project_folder_path, table):
    for dataset in ["galdieria", "micromonas", "phaeodactylum"]:
        for cell_id in range(1, 4):
            criteria = np.logical_and(table["species"].str.lower() == dataset, table["cell_id"] == cell_id)
            cut_table = table.loc[criteria, :].drop(columns=["species", "cell_id"])

            organelles_dir_path = os.path.join(project_folder_path, "data", dataset,
                                               "tables", f"cell{cell_id}-organelles")
            default_table_path = os.path.join(organelles_dir_path, "default.tsv")
            morphology_table_path = os.path.join(organelles_dir_path, "morphology.tsv")

            default_table = pd.read_csv(default_table_path, sep="\t")
            merged_table = default_table.merge(cut_table, how="left", on=["organelle"], validate="one_to_one")

            merged_table = merged_table[["label_id", "volume_micron_cubed",
                                         "surface_area_micron_squared", "ratio_percent"]]
            merged_table.to_csv(morphology_table_path, sep="\t", index=False)


def write_morphology_tables(project_folder_path, table_path):
    table = pd.read_csv(table_path, sep="\t")
    write_merged_dataset_tables(project_folder_path, table)
    write_unmerged_dataset_tables(project_folder_path, table)


if __name__ == '__main__':
    # project_folder = "C:\\Users\\meechan\\Documents\\Repos\\plankton-fibsem-project"
    # full_table = "C:\\Users\\meechan\\Documents\\collated_plankton_data.tsv"
    write_morphology_tables('.', './misc/collated_plankton_data.tsv')
