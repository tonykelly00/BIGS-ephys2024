
import pandas as pd
import numpy as np
import scanpy as sc

def generate_SelectMarkers(gene_data, metadata):
    
    # Assuming gene_data is your DataFrame with genes as rows and cells as columns
    adata = sc.AnnData(gene_data.T)

    # Ensure indices match
    metadata.index = adata.obs_names

    ## Assuming 'cell_type' is a column in adata.obs that contains cell type annotations
    adata.obs['cell_type'] = metadata['T-type Label'].astype('category')  

    ## Filter out groups with fewer than a minimum number of samples (e.g., 3)
    min_samples = 3
    cell_type_counts = adata.obs['cell_type'].value_counts()
    valid_cell_types = cell_type_counts[cell_type_counts >= min_samples].index
    adata = adata[adata.obs['cell_type'].isin(valid_cell_types)]

    # Log-transform the data
    sc.pp.log1p(adata)

    # # Find marker genes
    sc.tl.rank_genes_groups(adata, 'cell_type', method='t-test')

    # # Extract marker genes
    marker_genes = pd.DataFrame(adata.uns['rank_genes_groups']['names']).head(50)

    # # Save to CSV
    marker_genes.to_csv('../downloads/select_markers.csv', index=False)


def generate_gene_umap(gene_data):
    #method1
    # Read the CSV file
    marker_genes_for_umap = pd.read_csv("../downloads/select_markers.csv")

    # Melt the DataFrame to have a single column of gene names
    marker_genes_for_umap = marker_genes_for_umap.melt(value_name="Gene").dropna()["Gene"]

    # Ensure the gene names are unique
    marker_genes_for_umap = marker_genes_for_umap.drop_duplicates()

    # Check if the gene names exist in gene_data
    missing_genes = marker_genes_for_umap[~marker_genes_for_umap.isin(gene_data.index)]
    if not missing_genes.empty:
        print(f"Warning: The following genes are missing in gene_data: {missing_genes.tolist()}")

    return marker_genes_for_umap
