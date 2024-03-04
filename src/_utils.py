

import pandas as pd
import scanpy as sc
import numpy as np
from matplotlib import pyplot as plt
from anndata import AnnData
import sys

sys.path.insert(0, '/Rnd/utils/singlecell/marker_dict/')
from dotplot_profiling import modify_marker_list


#Vizgen function

def make_AnnData(cell_by_gene_path, meta_cell_path):
    
    cell_by_gene = pd.read_csv(cell_by_gene_path, index_col=0)
    meta_cell = pd.read_csv(meta_cell_path, index_col=0)

    meta_cell['barcodeCount'] = cell_by_gene.sum(axis=1)

    # initialize meta_gene
    meta_gene = pd.DataFrame(index=cell_by_gene.columns.tolist())

    # Align the cell id of cell_metadata to cell_by_gene
    cell_id = cell_by_gene.index.tolist()
    meta_cell = meta_cell.loc[cell_id]

    # Check again
    if (cell_by_gene.index == meta_cell.index).sum() == len(cell_by_gene.index):
        print('The indices in cell_by_gene and cell_metadata match.')
    else:
        print('The indices in cell_by_gene and cell_metadata do not match.')

    coordinates =np.stack((meta_cell['center_x'], meta_cell['center_y']), axis=1)

    ad = sc.AnnData(X=cell_by_gene.values, obs=meta_cell, var=meta_gene, obsm={"spatial": coordinates})
    return ad

# Generates dotplots with multiple dict of markers for cell type annotation

def dotplot_annotate(adata,marker_dict,list_dict,leiden):
    modify_marker_list(marker_dict, adata)
    for key, values in marker_dict.items():
        if key in list_dict:
            print(key)
            dp = sc.pl.dotplot(adata, marker_dict[key][0], groupby=leiden, return_fig=True,cmap ='Reds')
            dp.add_totals().show()  
    return



def plot_cell_type(adata,annotation):
    unique_categories = adata.obs[annotation].unique()

    # Determine number of rows and columns for subplots based on the number of unique categories
    nrows = int(len(unique_categories) ** 0.5)
    ncols = int(len(unique_categories) / nrows) + (len(unique_categories) % nrows > 0)

    fig, axes = plt.subplots(nrows=nrows, ncols=ncols, figsize=(20, 20))


    for i, category in enumerate(unique_categories):
        ax = axes[i // ncols, i % ncols]
        # Subset the data
        subset = adata[adata.obs[annotation] == category]
        # Plot

        
        sc.pl.spatial(adata[adata.obs['major_cell_type'].str.contains('Epithel')],add_outline =True, legend_loc=None, spot_size=30,show=False,palette='white',ax=ax,alpha=0)    
        sc.pl.spatial(subset,color=annotation ,legend_loc=None, spot_size=50,show=False, ax=ax, img_key="hires")
        ax.set_title(category, color="black")
        # Rotate the plot by 180 degrees
        #ax.invert_xaxis()
        ax.invert_yaxis()
    # If there are any remaining axes spots (in cases where the number of unique categories is not a perfect square), hide them
    for j in range(i+1, nrows*ncols):
        axes[j // ncols, j % ncols].axis('off')
        
        
        
def plot_spatial_gene(
    adata :AnnData,
    gene : str
    ) -> None:
    
    fig, ax = plt.subplots(figsize=(5, 5))  

    mask = adata[:, gene].X > 1 # Remove null expressiom
    adata_filtered = adata[mask]
    sc.pl.spatial(adata, legend_loc=None, spot_size=10,show=False, img_key="hires",palette='lightgray',ax=ax,alpha=0.5)    
    sc.pl.spatial(adata[adata.obs['major_cell_type'].str.contains('Epithel')],add_outline =True, legend_loc=None, spot_size=10,show=False, img_key="hires",palette='lightgray',ax=ax,alpha=0.5)    
    sc.pl.spatial(adata_filtered, color=[gene],legend_loc=None, spot_size=20,show=False, img_key="hires",cmap='gist_heat_r', ax=ax)
    plt.show()