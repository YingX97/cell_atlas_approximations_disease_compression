import os
import pathlib
import gc
import anndata
import cellxgene_census
import scquill


def compress_dataset_chunked(
    data_folder,
    dataset_id,
    include_neighborhood=False,
    overwrite=False,
    cache_adata=False,
    groupby=("tissue_general", "disease", "development_stage_general", "sex"),
    chunkby=("tissue_general", "cell_type"),
    optional_chunkby=("disease",),
):
    """Compress a dataset in chunks.

    The chunking is somewhat dynamic, based on the number of cells in the dataset. The logic is to chunk by
    avoid giant chunks at all costs since they are the ones bottlenecking RAM usage.
    """

    output_fn = f"{data_folder}/approximations/{dataset_id}.h5"
    if not overwrite and pathlib.Path(output_fn).exists():
        print(f"Already completed, skipping: {dataset_id}", flush=True)
        return dataset_id

    # The strategy is:
    # 1. Get the tissue metadata to split
    # 2. For each tissue, compress the data
    # 3. Merge the tissues back together
    # 4. Store to disk
    # 5. Add bells and whistles
    results_tissues = []

    print(f"Access CENSUS database to retrieve chunked data for dataset: {dataset_id}")
    # The most recent weekly release can be opened by the APIs by specifying census_version = "latest".
    with cellxgene_census.open_soma(census_version="latest") as census:

        obs = cellxgene_census.get_obs(
            census=census,
            organism="Homo sapiens",
            value_filter=f"dataset_id == '{dataset_id}'",
            column_names=[
                "cell_type",
                "tissue_general",
                "disease",
            ],
        )
        print(f"Access completed (tissue metadata): {dataset_id}")
        obs["c"] = 1

        # Set up the chunking
        chunkby_this = list(chunkby)
        chunk_cell_counts = obs.groupby(list(chunkby), observed=True).size()
        chunk_cell_counts = chunk_cell_counts[chunk_cell_counts >= 3]
        # If the chunking is too coarse, subchunk
        if chunk_cell_counts.max() > 5e4:
            chunkby_this = list(chunkby) + list(optional_chunkby)
            chunk_cell_counts = obs.groupby(list(chunkby_this), observed=True).size()
            chunk_cell_counts = chunk_cell_counts[chunk_cell_counts >= 3]

        chunks = chunk_cell_counts.index
        for chunk in chunks:

            print("Check adata cache")
            cache_fn = dataset_id
            for key, value in zip(chunkby_this, chunk):
                cache_fn += f"::{value}"
            cache_fn += ".h5ad"
            cache_fn = pathlib.Path(data_folder) / "__adata_cache__" / cache_fn
            if cache_fn.exists():
                print("Cache hit, reading h5ad file...")
                adata = anndata.read_h5ad(cache_fn)
            else:
                print("Cache miss, requesting adata from census...")
                obs_filter_string = f"(dataset_id == '{dataset_id}')"
                for key, value in zip(chunkby_this, chunk):
                    obs_filter_string += f" & ({key} == '{value}')"

                adata = cellxgene_census.get_anndata(
                    census=census,
                    organism="Homo sapiens",
                    obs_value_filter=obs_filter_string,
                    obs_column_names=[
                        "cell_type",
                        "assay",
                        "tissue",
                        "tissue_general",
                        "suspension_type",
                        "disease",
                        "sex",
                        "development_stage",
                    ],
                )
                if cache_adata:
                    print("Store adata cache file")
                    os.makedirs(cache_fn.parent, exist_ok=True)
                    adata.write(cache_fn)

            ncells = chunk_cell_counts[chunk]
            log_string = dataset_id
            for key, value in zip(chunkby_this, chunk):
                log_string += f", {value}"
            log_string += f" {ncells} cells"
            print(f"Access completed (counts): {log_string}")

            print("Postprocess cellxgene obs metadata")
            _postprocess_cellxgene_metadata(adata)

            print("Setting gene names as var_names")
            adata.var.set_index("feature_name", inplace=True, drop=False)
            # TODO: check that they are unique. They were unique in the one test I ran (Fabio)

            try:
                print(f"Build approximation for dataset: {dataset_id}")
                q = scquill.Compressor(
                    adata=adata,
                    celltype_column="cell_type",
                    additional_groupby_columns=groupby,
                    include_neighborhood=include_neighborhood,
                )
                q.prepare()
                q.compress()
                # FIXME: This misses the neighbourhoods!
                results_tissues.append(q.to_anndata())
                print(f"Successfully approximated: {log_string}")

                # Manual garbage management to reduce memory usage
                del q
                gc.collect()

            except Exception as e:
                print(
                    f"Error during compression for {log_string}: {e}",
                    flush=True,
                )
                raise

    print(f"Merging chunks from dataset: {dataset_id}")
    if len(chunks) == 1:
        result_dataset = results_tissues[0]
    else:
        result_dataset = anndata.concat(results_tissues)
    del results_tissues

    q = scquill.Compressor.from_anndata(
        adata=result_dataset,
        output_filename=output_fn,
    )
    q.store()

    del q, result_dataset
    gc.collect()

    print(f"Successfully approximated: {dataset_id}")

    return dataset_id


def _postprocess_cellxgene_metadata(adata: anndata.AnnData):
    adata.obs["development_stage_general"] = adata.obs["development_stage"].astype(
        object
    )
    stages = adata.obs["development_stage"].value_counts()
    stages = stages[stages >= 0].index
    for stage in stages:
        new_stage = _convert_development_stage(stage)
        adata.obs.loc[
            adata.obs["development_stage"] == stage, "development_stage_general"
        ] = new_stage


def _convert_development_stage(stage):
    if "mature" in stage.lower():
        return "adult"
    if "human adult" in stage.lower():
        return "adult"
    if "-year-old human" in stage.lower():
        return "adult"
    if "-month-old human" in stage.lower():
        return "child"
    if "fetal" in stage.lower():
        return "fetal"
    if "embryonic" in stage.lower():
        return "fetal"
