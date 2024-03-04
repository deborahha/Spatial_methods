import numpy as np
import scanpy as sc
import squidpy as sq
sc.settings.verbosity = 3   




def tumor_border(adata,radius):
    adata.obs.reset_index(inplace = True)
    sq.gr.spatial_neighbors(adata, coord_type="generic", radius=radius)
    Epithel_indices = set(np.where(adata.obs['major_cell_type'].str.contains('Epithel'))[0])
    Stroma_indices = set(np.where(adata.obs['major_cell_type'].isin(['Endothel','CAF']))[0])

    adata.obs['tumor_border'] = 'no'

    for idx in list(Epithel_indices):
        # Find neighbors of the current immune cell
        neighbors = adata.obsp['spatial_connectivities'][idx].nonzero()[1]

        if any(neighbor in Stroma_indices for neighbor in neighbors):
            adata.obs.loc[idx, 'tumor_border'] = 'yes'
    return adata