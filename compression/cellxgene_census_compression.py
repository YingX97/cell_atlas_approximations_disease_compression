import os
import scquill
import pathlib
import argparse
import numpy as np
import pandas as pd
import cellxgene_census
import multiprocess as mp
from datetime import datetime
from utils.guess_normalisation import (
    guess_unit_and_log,
    guess_unit_and_log_relaxed_method,
)
from utils.write_to_file import add_metadata, add_unit_and_log_to_compressed_dataset


def func(dataset_id):

    output_fn = f"data/approximations/{dataset_id}.h5"
    if pathlib.Path(output_fn).exists():
        print(f"Already completed, skipping: {dataset_id}", flush=True)
        return dataset_id

    print(f"Access CENSUS database to retrieve data for dataset: {dataset_id}")

    # census_version option: stable or latest.
    # The most recent weekly release can be opened by the APIs by specifying census_version = "latest".
    with cellxgene_census.open_soma(census_version="latest") as census:
        adata = cellxgene_census.get_anndata(
            census=census,
            organism="Homo sapiens",
            obs_value_filter=f"dataset_id == '{dataset_id}'",
            obs_column_names=[
                "assay",
                "cell_type",
                "tissue",
                "tissue_general",
                "suspension_type",
                "disease",
                "sex",
                "development_stage",
            ],
        )
    print(f"Access completed: {dataset_id}")

    print(f"Dataset size: {adata.shape}", flush=True)

    try:
        unit, log_transformed = guess_unit_and_log(adata)
        if unit is None:
            unit, log_transformed = guess_unit_and_log_relaxed_method(adata)

        if unit is None:
            exit(f"Cannot guess normalisation for {dataset_id}")

        normalisation = unit if not log_transformed else unit + "+logp1"

        print(f"Build approximation for dataset: {dataset_id}")
        q = scquill.Compressor(
            adata=adata,
            celltype_column="cell_type",
            output_filename=output_fn,
            additional_groupby_columns=(
                "tissue",
                "tissue_general",
                "disease",
                "sex",
                "development_stage",
            ),
            configuration={"normalisation": normalisation} if normalisation else {},
        )
        q()

        print(f"Successfully approximated: {dataset_id}")

        # first add the unit and whether it is logged based on the uncompressed adata object
        add_unit_and_log_to_compressed_dataset(output_fn, unit, log_transformed)

        # second, use the compressed data, and summarise the metadata and add to the file's attribute
        add_metadata(output_fn)

    except Exception as e:
        print(f"Error during compression for dataset {dataset_id}: {e}", flush=True)
        raise

    if pathlib.Path(output_fn).exists():
        print(f"Successfully created: {output_fn}", flush=True)
    else:
        print(f"Failed to create: {output_fn}", flush=True)
        raise RuntimeError(f"Failed to create output file for dataset {dataset_id}")

    return dataset_id


def pool_callback(dataset_id):
    print(f"Finished pool job for dataset: {dataset_id}", flush=True)
    # HERE COMES THE MERGIN OF THE TISSUES


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--threads",
        type=int,
        default=1,
        help="Number of processes to use in parallel",
    )
    parser.add_argument(
        "--test",
        action="store_true",
        help="Test only the first dataset, on 1 core",
    )
    parser.add_argument(
        "--dry",
        action="store_true",
        help="Dry run, don't call anything except the dataset_id generation",
    )
    args = parser.parse_args()

    # Only compress data with normal baseline
    # multi_condition_datasets = pd.read_csv(
    #     "../static/multi_condition_datasets.csv", index_col=0
    # )

    today = datetime.today().strftime("%d_%b_%Y")
    datasets_fn = pathlib.Path(f"data/cellxgene_census_dataset_id_{today}.csv")

    if not pathlib.Path("data").exists():
        os.mkdir("data")
    if not pathlib.Path("data/approximations").exists():
        os.mkdir("data/approximations")

    if not datasets_fn.exists():
        print("Access CENSUS database to retrieve dataset ids")
        with cellxgene_census.open_soma() as census:
            dataset_ids = (
                census["census_data"]["homo_sapiens"]
                .obs.read(column_names=["dataset_id"])
                .concat()
            )
            # They are already sorted by size, now invert it
            dataset_id_count = (
                dataset_ids["dataset_id"].to_pandas().value_counts()[::-1]
            )
        print("Access completed")

        dataset_id_count.to_csv(datasets_fn, sep=",")
        print("Done")

    print("Load dataset ids")
    dataset_id_count = pd.read_csv(datasets_fn, index_col=0)["count"]
    dataset_ids = dataset_id_count.index.values

    if args.dry:
        print(f"Dry run with {args.threads} processes and test={args.test}")
    elif args.test:
        for dataset_id in dataset_ids:
            func(dataset_id)
            break
    elif args.threads == 1:
        ndata = len(dataset_ids)
        for i, dataset_id in enumerate(dataset_ids):
            print(f"Start dataset {i+1} / {ndata}: {dataset_id}", flush=True)
            func(dataset_id)

    else:
        with mp.Pool(args.threads) as pool:
            # Get promise for good coding style
            async_res = pool.map_async(func, dataset_ids, callback=pool_callback)
            # Wait for the end on this process
            async_res.wait()
