from datetime import datetime
import cellxgene_census
import h5py
import json
import numpy as np
import scquill
import hashlib
import uuid

# load the cellxgene census dataset information
# this will include collection_doi, collection_name etc
# add these metadata to a file during compression
with cellxgene_census.open_soma(census_version='latest') as census:
    CENSUS_DATASET_INFO = census["census_info"]["datasets"].read().concat().to_pandas()
    CENSUS_DATASET_INFO = CENSUS_DATASET_INFO.set_index("dataset_id").to_dict('index')


def add_to_manifest(dataset_id, information):
    # Implementation of add_to_manifest (same as provided)
    today = datetime.today().strftime('%d_%b_%Y')
    try:
        with open(f"data/compression_manifest_{today}.json") as f:
            manifest = json.load(f)
    except:
        manifest = {}
    
    metadata = manifest[dataset_id] if dataset_id in manifest else CENSUS_DATASET_INFO[dataset_id]
    manifest[dataset_id] = {**metadata, **information}
    with open(f"data/compression_manifest_{today}.json", 'w') as f:
        json.dump(manifest, f, indent=2)
        


def add_unit_and_log_to_compressed_dataset(compressed_file_path, unit, log_transformed):
    # Implementation of add_unit_and_log_to_compressed_dataset (same as provided)
    dataset_id = compressed_file_path.split("/")[-1].replace(".h5", "")
    add_to_manifest(dataset_id, {"unit": unit, "log_transformed": log_transformed})
    with h5py.File(compressed_file_path, "a") as compressed_h5:
        compressed_h5.attrs["unit"] = unit
        compressed_h5.attrs["log_transformed"] = log_transformed


def add_metadata(compressed_file_path):
    # Implementation of add_metadata (same as provided)
    dataset_id = compressed_file_path.split("/")[-1].replace(".h5", "")
    app = scquill.Approximation()
    app = app.read_h5(compressed_file_path)
    adata = app.to_anndata(
        groupby=(
            "cell_type",
            "tissue",
            "tissue_general",
            "disease",
            "sex",
            "development_stage",
        )
    )
    compressed_obs = adata.obs
    observations = {}
    for combination in compressed_obs.index:
        cell_type = combination.split("\t")[0]
        tissue = combination.split("\t")[1]
        tissue_general = combination.split("\t")[2]
        disease = combination.split("\t")[3]
        sex = combination.split("\t")[4]
        development_stage = combination.split("\t")[5]
        if disease in observations:
            observations[disease]["cell_type"].append(cell_type)
            observations[disease]["tissue"].append(tissue)
            observations[disease]["tissue_general"].append(tissue_general)
            observations[disease]["sex"].append(sex)
            observations[disease]["development_stage"].append(development_stage)
        else:
            observations[disease] = {}
            observations[disease]["cell_type"] = [cell_type]
            observations[disease]["tissue"] = [tissue]
            observations[disease]["tissue_general"] = [tissue_general]
            observations[disease]["sex"] = [sex]
            observations[disease]["development_stage"] = [development_stage]
    has_normal_baseline = "normal" in observations
    cell_type = []
    sex = []
    tissue = []
    tissue_general = []
    development_stage = []
    unique_ids = []
    for disease in observations:
        if disease != "normal":
            cell_type += observations[disease]["cell_type"]
            tissue += observations[disease]["tissue"]
            tissue_general += observations[disease]["tissue_general"]
            sex += observations[disease]["sex"]
            development_stage += observations[disease]["development_stage"]
            item = {
                "dataset_id": dataset_id,
                "disease": disease,
                "cell_type": observations[disease]["cell_type"],
                "tissue": observations[disease]["tissue"],
                "tissue_general": observations[disease]["tissue_general"],
                "sex": observations[disease]["sex"],
                "development_stage": observations[disease]["development_stage"],
            }
            stringified = json.dumps(item)
            m = hashlib.md5()
            m.update(stringified.encode("utf-8"))
            unique_ids.append(str(uuid.UUID(m.hexdigest())))
    cell_type = list(set(cell_type))
    tissue = list(set(tissue))
    tissue_general = list(set(tissue_general))
    development_stage = list(set(development_stage))
    sex = list(set(sex))

    add_to_manifest(
        dataset_id,
        {
            "ids": unique_ids,
            "cell_type": cell_type,
            "tissue_general": tissue_general,
            "disease": [disease for disease in observations if disease != "normal"],
            "sex": sex,
            "development_stage": development_stage,
            "has_normal_baseline": has_normal_baseline,
        },
    )
    with h5py.File(compressed_file_path, "a") as compressed_h5:
        compressed_h5.attrs["ids"] = np.array(
            [uid.encode("UTF-8", "ignore") for uid in unique_ids]
        )
        compressed_h5.attrs["dataset_id"] = dataset_id
        compressed_h5.attrs["cell_type"] = np.array(
            [c.encode("UTF-8", "ignore") for c in cell_type]
        )
        compressed_h5.attrs["tissue_general"] = np.array(
            [t.encode("UTF-8", "ignore") for t in tissue_general]
        )
        compressed_h5.attrs["disease"] = np.array(
            [
                disease.encode("UTF-8", "ignore")
                for disease in observations
                if disease != "normal"
            ]
        )
        compressed_h5.attrs["sex"] = np.array(
            [s.encode("UTF-8", "ignore") for s in sex]
        )
        compressed_h5.attrs["development_stage"] = np.array(
            [d.encode("UTF-8", "ignore") for d in development_stage]
        )
        compressed_h5.attrs["has_normal_baseline"] = has_normal_baseline
