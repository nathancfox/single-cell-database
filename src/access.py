import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import os
import h5py as h5
import loompy as lp
import scanpy as sc
import pandas as pd
import global_constants as GC
import external_metadata as em__

def get_loom_filename(uuid):
    df = em__.get_as_dataframe()
    if uuid not in list(df['uuid']):
        raise ValueError('uuid is not valid!')
    filename = df[df['uuid'] == uuid]['file_location'].iloc[0]
    if not os.path.exists(filename):
        raise AssertionError('Retrieved filename does not exist!')
    filename = os.path.join(filename, 'expr_mat.loom')
    return(filename)

def get_h5_conn(uuid):
    lfile = h5.File(get_loom_filename(uuid), 'r')
    return(lfile)

def get_loom_conn(uuid):
    lfile = lp.connect(get_loom_filename(uuid), 'r')
    return(lfile)

def get_anndata(uuid, **kwargs):
    # Handles the loom convention that genes may be named
    # in the 'Gene' or 'Accession' row attribute.
    if 'var_names' not in kwargs.keys():
        lfile = get_h5_conn(uuid)
        if 'Gene' in lfile['row_attrs'].keys():
            pass
        elif 'Accession' in lfile['row_attrs'].keys():
            kwargs['var_names'] = 'Accession'
        else:
            pass
    adata = sc.read_loom(get_loom_filename(uuid), **kwargs)
    return(adata)

def main():
    pass

if __name__ == '__main__':
    main()