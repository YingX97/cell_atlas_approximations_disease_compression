import pathlib
import gc
import anndata
import cellxgene_census
import scquill

from utils.guess_normalisation import (
    guess_unit_and_log,
    guess_unit_and_log_relaxed_method,
)
from utils.write_to_file import (
    add_metadata,
    add_unit_and_log_to_compressed_dataset,
)


def compress_dataset(
    data_folder,
    dataset_id,
    include_neighborhood=False,
    overwrite=False,
):
    output_fn = f"{data_folder}/approximations/{dataset_id}.h5"
    if not overwrite and pathlib.Path(output_fn).exists():
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

    print("Setting gene names as var_names")
    adata.var.set_index("feature_name", inplace=True, drop=False)
    # TODO: check that they are unique. They were unique in the one test I ran (Fabio)

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
            include_neighborhood=include_neighborhood,
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


def compress_dataset_chunked(
    data_folder,
    dataset_id,
    include_neighborhood=False,
    overwrite=False,
):
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
                "tissue_general",
            ],
        )
        print(f"Access completed (tissue metadata): {dataset_id}")

        # Extract the tissues
        tissue_counts = obs["tissue_general"].value_counts()
        tissues = tissue_counts[tissue_counts > 0].index

        for tissue in tissues:
            adata = cellxgene_census.get_anndata(
                census=census,
                organism="Homo sapiens",
                obs_value_filter=f"(dataset_id == '{dataset_id}') & (tissue_general == '{tissue}')",
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
            print(
                f"Access completed (counts): {dataset_id}, {tissue}, {adata.n_obs} cells"
            )

            if adata.n_obs < 50:
                continue

            print("Setting gene names as var_names")
            adata.var.set_index("feature_name", inplace=True, drop=False)
            # TODO: check that they are unique. They were unique in the one test I ran (Fabio)

            try:
                unit, log_transformed = guess_unit_and_log(adata)
                if unit is None:
                    unit, log_transformed = guess_unit_and_log_relaxed_method(adata)

                if unit is None:
                    exit(f"Cannot guess normalisation for {dataset_id}, {tissue}")

                normalisation = unit if not log_transformed else unit + "+logp1"

                print(f"Build approximation for dataset: {dataset_id}, {tissue}")
                q = scquill.Compressor(
                    adata=adata,
                    celltype_column="cell_type",
                    additional_groupby_columns=(
                        "tissue",
                        "tissue_general",
                        "disease",
                        "sex",
                        "development_stage",
                    ),
                    configuration=(
                        {"normalisation": normalisation} if normalisation else {}
                    ),
                    include_neighborhood=include_neighborhood,
                )
                q.prepare()
                q.compress()
                # FIXME: This misses the neighbourhoods!
                results_tissues.append(q.to_anndata())
                print(f"Successfully approximated: {dataset_id}, {tissue}")

                # Manual garbage management to reduce memory usage
                del q
                gc.collect()

            except Exception as e:
                print(
                    f"Error during compression for dataset {dataset_id}: {e}",
                    flush=True,
                )
                raise

    print(f"Merging tissues from dataset: {dataset_id}")
    if len(tissues) == 1:
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

    # first add the unit and whether it is logged based on the uncompressed adata object
    add_unit_and_log_to_compressed_dataset(output_fn, unit, log_transformed)

    # second, use the compressed data, and summarise the metadata and add to the file's attribute
    add_metadata(output_fn)
    print(f"Successfully approximated: {dataset_id}")

    return dataset_id
