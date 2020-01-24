"""Abstracted access functions to the database.

A collection of functions that provide access to the
database with minimal knowledge of the inner workings.
In general, the external metadata can be gotten with
0 prior knowledge, which can be parsed for the desired
UUID(s). All the other access functions only require a
UUID.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import os
import h5py as h5
import loompy as lp
import scanpy as sc
import numpy as np
import pandas as pd
from . import external_metadata as em__
from . import internal_metadata as im__
from . import global_constants as GC

def get_loom_filename(uuid):
    """Get loom filename for dataset entry.

    Given a UUID, get the filepath to the associated
    loom file from the external metadata file. If the
    UUID is not in the external metadata file, the
    corresponding folder is searched for in the database.
    This allows this function to be used during a new
    entry, when the loom file is created, but the entry
    hasn't been entered into the external metadata yet.

    Args:
        uuid: String. UUID of the desired dataset
    
    Returns:
        String containing the filepath of the loom file
        associated with that UUID.

    Raises:
        ValueError: Raised if uuid is not in the external
            metadata file and is also not a folder with an
            "expr_mat.loom" file in the database.
        AssertionError: Raised if the filepath to be returned
            does not exist.
    """
    df = em__.get_as_dataframe()
    if uuid not in list(df['uuid']):
        for dir_entry in os.scandir(GC._PATH_TO_DATABASE):
            if (dir_entry.name == uuid
                    and dir_entry.is_dir()
                    and 'expr_mat.loom' in os.listdir(dir_entry.path)):
                return(os.path.join(dir_entry.path, 'expr_mat.loom'))
        raise ValueError('uuid is not valid!')
    else:
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
        with get_h5_conn(uuid) as lfile:
            if 'Accession' in lfile['row_attrs'].keys():
                kwargs['var_names'] = 'Accession'
            elif 'Gene' in lfile['row_attrs'].keys():
                kwargs['var_names'] = 'Gene'
            else:
                # By passing a column that doesn't exist, the constructor
                # silently creates a numbered index
                kwargs['var_names'] = 'PLACEHOLDER_FOR_NUMBERED_INDEX'
    adata = sc.read_loom(get_loom_filename(uuid), **kwargs)
    # Accession and Gene are handled separately because they
    # aren't handled with the global constants, and so they
    # can't be sorted by their index.
    add_acc = False
    add_gene = False
    if 'Accession' in list(adata.var.columns):
        add_acc = True
    if 'Gene' in list(adata.var.columns):
        add_gene = True
    var_columns = list(filter(lambda x: x not in ['Accession', 'Gene'],
                              adata.var.columns))
    var_column_order = sorted(var_columns,
                              key = lambda x: GC._IMU_GENE_COLUMN_INDEX[x])
    # Prepended in reverse order to ensure that
    # Accession will come before Gene
    if add_gene:
        var_column_order = ['Gene'] + var_column_order
    if add_acc:
        var_column_order = ['Accession'] + var_column_order
    adata.var = adata.var[var_column_order]
    adata.obs = adata.obs[sorted(adata.obs.columns,
                                 key = lambda x: GC._IMU_CELL_COLUMN_INDEX[x])]
    adata.uns['batch_key'] = get_batch_key(uuid)
    adata.uns['cell_author_annot'] = get_cell_author_annot(uuid)
    adata.uns['gene_author_annot'] = get_gene_author_annot(uuid)
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
    """Get cell IDs from a dataset.

    Get the cell IDs from a specific dataset.

    Args:
        uuid: String. The UUID of the desired dataset.

    Returns:
        A 1-D numpy array of strings holding the cell IDs.

    Raises: None
    """
    with get_h5_conn(uuid) as lfile:
        cell_ids = np.array(lfile['col_attrs/CellID'])
        return(cell_ids)

def get_gene_ids(uuid, accession = False):
    """Get gene IDs from a dataset.

    Args:
        uuid: String. The UUID of the desired dataset.
        accession: Boolean. If True, the "Accession" field
            will be returned. If False, the "Gene" field
            will be returned. If the chosen field is not
            available and the other is, a warning will be
            printed and the other will be returned.
        
    Returns:
        A 1-D numpy array of strings holding the gene IDs.

    Raises:
        AssertionError: If the dataset does not have either a
            "Gene" field or an "Accession" field.
    """
    with get_h5_conn(uuid) as lfile:
        if accession:
            if 'Accession' in lfile['row_attrs'].keys():
                gene_ids = np.array(lfile['row_attrs/Accession']) 
            elif 'Gene' in lfile['row_attrs'].keys():
                print('Warning! "Accession" not available. Returning \"Gene\" '
                    'instead.')
                gene_ids = np.array(lfile['row_attrs/Gene']) 
            else:
                raise AssertionError(f'Dataset {uuid} does not have a \"Gene\" or '
                                    f'an \"Accession\" row attribute!')
        else:
            if 'Gene' in lfile['row_attrs'].keys():
                gene_ids = np.array(lfile['row_attrs/Gene']) 
            elif 'Accession' in lfile['row_attrs'].keys():
                print('Warning! "Gene" not available. Returning \"Accession\" '
                    'instead.')
                gene_ids = np.array(lfile['row_attrs/Accession']) 
            else:
                raise AssertionError(f'Dataset {uuid} does not have a \"Gene\" or '
                                    f'an \"Accession\" row attribute!')
        return(gene_ids)

def get_column_description(uuid, column, var = 'cell', metadata = 'universal'):
    """Get the description for an internal metadata column.

    This is a wrapper around methods from internal_metadata.py
    Get the HDF5 attribute "description" from a column in
    the internal universal or author_annot metadata for a
    given dataset.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired column.
        var: String. Must be "cell" or "gene". Indicates
            the metadata to be queried.
        metadata: String. Must be "universal" or "author_annot".
            Indicates the metadata to be queried.

    Returns:
        Whatever is at the "description" HDF5 attribute of the
        requested column. If no description is available, a
        warning will be printed and the method will return None.

    Raises:
        ValueError: If var or metadata is passed an invalid value.
    """
    if metadata == 'universal':
        if var == 'cell':
            return(im__.get_cell_univ_col_desc(uuid, column))
        elif var == 'gene':
            return(im__.get_gene_univ_col_desc(uuid, column))
        else:
            raise ValueError('var must be \"cell\" or \"gene\"!')
    elif metadata == 'author_annot':
        if var == 'cell':
            return(im__.get_cell_aa_col_desc(uuid, column))
        elif var == 'gene':
            return(im__.get_gene_aa_col_desc(uuid, column))
        else:
            raise ValueError('var must be \"cell\" or \"gene\"!')
    else:
        raise ValueError('metadata must be \"universal\" or \"author_annot\"!')

def get_expr_mat(uuid, matrix = 'matrix'):
    """Get an expression matrix as a numpy array.

    From the given dataset, extract an expression matrix from it
    as a numpy array.

    Args:
        uuid: String. The UUID of the desired dataset.
        matrix: String. The name of the desired expression matrix.
            Must be 'matrix' or the name of an expression matrix
            under lfile['layers/'].
        
    Returns:
        A numpy array containing the expression matrix.
        Genes are rows, and cells are columns.

    Raises:
        ValueError: If the requested matrix does not exist.
    """
    with get_h5_conn(uuid) as lfile:
        if matrix == 'matrix':
            mat = np.array(lfile['matrix'])
        else:
            if matrix not in lfile['layers'].keys():
                raise ValueError('matrix must be \"matrix\" or a '
                                 'valid layer name!')
            else:
                mat = np.array(lfile[f'layers/{matrix}'])
        return(mat)

def get_batch_key(uuid):
    """Get the batch_key from a dataset.

    This is a wrapper around internal_metadata.get_cell_batch_key()

    Args:
        uuid: String. The UUID of the desired dataset.

    Returns:
        The "batch_key" string stored in the dataset.
    
    Raises:
        AssertionError: If the "batch" column or the
            "batch_key" attribute doesn't exist.
    """
    batch_key = im__.get_cell_batch_key(uuid)
    return(batch_key)

def get_expr_mat_names(uuid):
    """Get the expression matrix names from a dataset.

    Get the names of all expression matrices in a dataset.
    The first entry will always be "matrix", referring to
    the expression matrix stored under lfile['matrix'],
    instead of under lfile['layers'].

    Args:
        uuid: String. The UUID Of the desired dataset.

    Returns:
        A 1-D numpy array of strings containing
        the names of the matrices.

    Raises: None:
    """
    with get_h5_conn(uuid) as lfile:
        names = ['matrix']
        for k in lfile['layers'].keys():
            names.append(k)
        names = np.array(names)
        return(names)

def main():
    pass

if __name__ == '__main__':
    main()