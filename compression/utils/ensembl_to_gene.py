import pandas as pd
import requests

def fetch_gene_pairs(species):
    """
    Fetches a list of gene pairs (Ensembl Gene ID and external gene name) for a given species from the Ensembl BioMart service.

    Parameters:
    species (str): The species identifier (e.g., "hsapiens" for human, "mmusculus" for mouse).
    
    Returns:
    pd.DataFrame: A pandas DataFrame containing two columns:
        - 'ensembl_id': Ensembl Gene ID
        - 'gene_id': External (common) gene name
    """
    
    server = "http://www.ensembl.org/biomart/martservice"
    query = f"""<?xml version="1.0" encoding="UTF-8"?>
    <!DOCTYPE Query>
    <Query virtualSchemaName="default" formatter="TSV" header="0" uniqueRows="1" count="" datasetConfigVersion="0.6">
        <Dataset name="{species}_gene_ensembl" interface="default">
            <Attribute name="ensembl_gene_id" />
            <Attribute name="external_gene_name" />
        </Dataset>
    </Query>"""
    
    response = requests.post(server, data={'query': query})
    response.raise_for_status()
    
    # Convert the response text to a pandas dataframe
    data = [line.split('\t') for line in response.text.strip().split('\n')]
    df = pd.DataFrame(data, columns=['ensembl_id', 'gene_id'])
    
    return df

# Fetch data for human and mouse
human_df = fetch_gene_pairs("hsapiens")
mouse_df = fetch_gene_pairs("mmusculus")

# Save data to CSV files
human_df.to_csv("../../static/human_gene_pairs.csv", index=False)
mouse_df.to_csv("../../static/mouse_gene_pairs.csv", index=False)

print("Data saved successfully!")
