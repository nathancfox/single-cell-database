import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import os
import h5py as h5
import loompy as lp
import scanpy as sc
import numpy as np
import pandas as pd
import external_metadata as em__
import internal_metadata as im__

def get_loom_filename(uuid):
    """Get loom filename for dataset entry.

    Given a UUID, get the filepath to the associated
    loom file from the external metadata file.

    Args:
        uuid: String. UUID of the desired dataset
    
    Returns:
        String containing the filepath of the loom file
        associated with that UUID.

    Raises:
        ValueError: Raised if uuid is not in the external
            metadata file.
        AssertionError: Raised if the filepath to be returned
            does not exist.
    """
    df = em__.get_as_dataframe()
    if uuid not in list(df['uuid']):
        raise ValueError('uuid is not valid!')
    filename = df[df['uuid'] == uuid]['file_location'].iloc[0]
    if not os.path.exists(filename):
        raise AssertionError('Retrieved filename does not exist!')
    filename = os.path.join(filename, 'expr_mat.loom')
    return(filename)

def get_h5_conn(uuid, write = False):
    """Get h5py file connection for dataset loom file.

    Given a UUID, get an h5py File connection object
    for the associated loom file.

    Args:
        uuid: String. UUID of the desired dataset
        write: Boolean. If True, the File object
            will have write access. Note that this
            is a temporary feature during development
            and will not be available in the final version.

    Returns:
        h5py File object to the desired loom file

    Raises: None
    """
    if write:
        lfile = h5.File(get_loom_filename(uuid), 'r+')
    else:
        lfile = h5.File(get_loom_filename(uuid), 'r')
    return(lfile)

def get_loom_conn(uuid):
    """Get loompy file connection for dataset loom file.

    Given a UUID, get a loompy file connection object
    for the associated loom file.

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        loompy LoomConnection object to the desired loom file

    Raises: None
    """
    lfile = lp.connect(get_loom_filename(uuid), 'r')
    return(lfile)

def get_anndata(uuid, **kwargs):
    """Get an AnnData object from a dataset.

    Given a UUID, loads the associated loom file
    completely into memory into a scanpy AnnData object.

    Args:
        uuid: String. UUID of the desired dataset
        **kwargs: Keyword arguments passed to
            scanpy.read_loom()
    
    Returns:
        AnnData object holding the dataset from the
        desired loom file.
    
    Raises: None
    """
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
    lfile.close()
    adata = sc.read_loom(get_loom_filename(uuid), **kwargs)
    return(adata)

def get_cell_univ(uuid, keep_missing = True):
    """Get the cell universal metadata.

    Get the cell universal metadata from the given dataset
    as a Pandas DataFrame. Wrapper for the actual method
    in internal_metadata.py

    Args:
        uuid: String. UUID of the desired dataset
        keep_missing: bool. If True, missing columns will
            be retained as columns of all -1 or "-1".
            If False, they are dropped from the returned DataFrame.
    
    Returns:
        Pandas DataFrame with the same number of rows as cells
        in the dataset and where each column is a universal
        internal metadata field.

    Raises: None
    """
    data = im__.get_cell_int_md_univ(uuid, keep_missing)
    return(data)

def get_gene_univ(uuid, keep_missing = True):
    """Get the gene universal metadata.

    Get the gene universal metadata from the given dataset
    as a Pandas DataFrame. Wrapper for the actual method
    in internal_metadata.py

    Args:
        uuid: String. UUID of the desired dataset
        keep_missing: bool. If True, missing columns will
            be retained as columns of all -1 or "-1".
            If False, they are dropped from the returned DataFrame.
    
    Returns:
        Pandas DataFrame with the same number of rows as genes
        in the dataset and where each column is a universal
        internal metadata field.

    Raises: None
    """
    data = im__.get_gene_int_md_univ(uuid, keep_missing)
    return(data)

def get_cell_author_annot(uuid):
    """Get the cell author-annotated metadata.

    Get the cell author-annotated metadata from the given
    dataset as a Pandas DataFrame. Wrapper for the
    actual method in internal_metadata.py

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as cells
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    data = im__.get_cell_int_md_author_annot(uuid)
    return(data)

def get_gene_author_annot(uuid):
    """Get the gene author-annotated metadata.

    Get the gene author-annotated metadata from the given
    dataset as a Pandas DataFrame. Wrapper for the
    actual method in internal_metadata.py

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as genes
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    data = im__.get_gene_int_md_author_annot(uuid)
    return(data)

def get_extern_md():
    """Gets external metadata as pandas DataFrame."""
    return(em__.get_as_dataframe())

def uuid_to_row(uuid, columns = None):
    """Get values of a row, indexed by UUID.

    This is just a wrapper around
    external_metadata.uuid_to_row()

    Given a certain UUID, retrieves that row from
    the external metadata and returns the requested
    columns as a pandas Series. If columns is None,
    the entire row is returned. If columns is a list
    of column names, the values are returned in that
    order. Invalid columns do not raise an Exception,
    however, a warning is generated and all values
    for those columns will be 'NOT_A_VALID_COLUMN'.

    Args:
        uuid_key: String. UUID used to look up the
            row requested. Must be in the 32 hexadecimal
            digits format created by str(uuid), where
            uuid is a python uuid.UUID class.
            e.g. "12345678-1234-5678-1234-567812345678"
        columns: list of strings. Each string should be
            a column from the external metadata file
            to include in the returned values. If None,
            the entire row will be returned.
    
    Returns:
        pandas Series containing the values requested
        that correspond to the uuid_key. Index should
        be the columns variable, if passed. Otherwise,
        will be the columns of the dataframe.

    Raises:
        IndexError: if the uuid_key is not a valid lookup.
    """
    return(em__.uuid_to_row(uuid, columns))

def get_cell_ids(uuid):
    lfile = get_h5_conn(uuid)
    cell_ids = np.array(lfile['col_attrs/CellID'])
    lfile.close()
    return(cell_ids)

def get_gene_ids(uuid, accession = False):
    lfile = get_h5_conn(uuid)
    if accession:
        if 'Accession' in lfile['row_attrs'].keys():
            gene_ids = np.array(lfile['row_attrs/Accession']) 
        elif 'Gene' in lfile['row_attrs'].keys():
            print('Warning! Accession not available. Returning \"Gene\" '
                  'instead.')
            gene_ids = np.array(lfile['row_attrs/Gene']) 
        else:
            raise AssertionError(f'Dataset {uuid} does not have a \"Gene\" or '
                                 f'an \"Accession\" row attribute!')
    else:
        if 'Gene' in lfile['row_attrs'].keys():
            gene_ids = np.array(lfile['row_attrs/Gene']) 
        elif 'Accession' in lfile['row_attrs'].keys():
            print('Warning! Gene not available. Returning \"Accession\" '
                  'instead.')
            gene_ids = np.array(lfile['row_attrs/Accession']) 
        else:
            raise AssertionError(f'Dataset {uuid} does not have a \"Gene\" or '
                                 f'an \"Accession\" row attribute!')
    lfile.close()
    return(gene_ids)

def main():
    pass

if __name__ == '__main__':
    main()