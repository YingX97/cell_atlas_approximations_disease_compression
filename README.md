# Cell Atlas Approximations Disease Compression

This repository contains scripts for compressing datasets from the cellxgene census using the [`scquill`](https://github.com/fabilab/scquill/tree/main) package. 
## Repository Structure

- `compression/`: Contains the main compression scripts and utilities.
  - `cellxgene_census_compression.py`: Main script for compressing cellxgene census datasets.
  - `utils/`: Utility functions for the compression process.
    - `ensembl_to_gene.py`: Converts Ensembl IDs to gene names.
    - `guess_normalisation.py`: Determines the normalization method used in the dataset.
    - `write_to_file.py`: Handles writing metadata(dataset title, collections...etc) to files.

- `static/`: Contains static data files used in the compression process.
  - `human_gene_pairs.csv`: Gene pair information for human.
  - `mouse_gene_pairs.csv`: Gene pair information for mouse.
  - `multi_condition_datasets.csv`: Contains dataset IDs that have multiple experimental conditions.

- `requirements.txt`: Lists all Python dependencies for the project.

## Getting Started

These instructions will guide you through setting up your project environment and running the compression scripts.

### Prerequisites

- Python 3.11s
- pip (Python package installer)

### 1. Clone the Repository

First, clone this repository to your local machine using git:

```bash
git clone https://github.com/YingX97/cell_atlas_approximations_disease_compression
cd cell_atlas_approximations_disease_compression
```

### 2. Set up Python virtual environment

It is recommended to use a virtual environment to avoid conflicts with other projects or system packages:

```bash
python3.11 -m venv venv
source venv/bin/activate  # On Windows use `venv\Scripts\activate`
```

### 3. Install required packages

Install all required packages using pip:

```bash
pip install -r requirements.txt
```

### 4. Running the compression script

```bash
cd compression
python3 cellxgene_census_compression.py
```