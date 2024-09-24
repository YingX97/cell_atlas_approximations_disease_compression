import os
import pathlib
import argparse
import numpy as np
import pandas as pd
import cellxgene_census
import multiprocess as mp
import concurrent.futures
from datetime import datetime

from utils.dbs.cellxgene_census import (
    compress_dataset as _compress_dataset,
    compress_dataset_chunked as _compress_dataset_chunked,
)


data_folder = pathlib.Path("__file__").absolute().parent.parent / "data"


def compress_dataset(dataset_id, chunked=True, include_neighborhood=False):
    """Worker routine to compress the dataset"""
    if chunked:
        return _compress_dataset_chunked(
            data_folder, dataset_id, include_neighborhood=include_neighborhood
        )
    else:
        return _compress_dataset(
            data_folder, dataset_id, include_neighborhood=include_neighborhood
        )


def pool_callback(dataset_id):
    print(f"Finished pool job for dataset: {dataset_id}", flush=True)


if __name__ == "__main__":

    parser = argparse.ArgumentParser()
    parser.add_argument(
        "--max-datasets",
        type=int,
        default=-1,
        help="Maximum number of datasets to process",
    )
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
    parser.add_argument(
        "--multiprocess",
        action="store_true",
        help="Use multiprocessing instead of threading",
    )
    args = parser.parse_args()

    today = datetime.today().strftime("%d_%b_%Y")
    datasets_fn = pathlib.Path(f"{data_folder}/cellxgene_census_dataset_id_{today}.csv")

    if not pathlib.Path("data").exists():
        os.mkdir("data")
    if not pathlib.Path(f"{data_folder}/approximations").exists():
        os.mkdir(f"{data_folder}/approximations")

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

    if args.max_datasets != -1:
        print(f"Restricting to {args.max_datasets} datasets")
        dataset_ids = dataset_ids[: args.max_datasets]

    if args.dry:
        print(f"Dry run with {args.threads} processes and test={args.test}")
    elif args.test:
        for dataset_id in dataset_ids:
            compress_dataset(dataset_id)
            break
    elif args.threads == 1:
        ndata = len(dataset_ids)
        for i, dataset_id in enumerate(dataset_ids):
            print(f"Start dataset {i+1} / {ndata}: {dataset_id}", flush=True)
            compress_dataset(dataset_id)

    elif not args.multiprocess:
        with concurrent.futures.ThreadPoolExecutor(
            max_workers=args.threads
        ) as executor:
            futures = []
            for i, dataset_id in enumerate(dataset_ids):
                print(f"Start dataset {i+1} / {len(dataset_ids)}: {dataset_id}")
                future = executor.submit(compress_dataset, dataset_id)
                future.add_done_callback(pool_callback)
                futures.append(future)
            # Wait for the end on this process
            concurrent.futures.wait(futures)
    else:
        with mp.Pool(args.threads) as pool:
            # Get promise for good coding style
            async_res = pool.map_async(
                compress_dataset, dataset_ids, callback=pool_callback
            )
            # Wait for the end on this process
            async_res.wait()
