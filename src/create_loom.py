"""Functions to parse expression data into a new loom file.

The main function is create_loom_file(). This creates a new
loom file for a new entry from a standard set of arguments.
All other functions exist to create that standard set of
arguments from a variety of expression data formats.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import os
import loompy as lp
import numpy as np
import pandas as pd
import scipy
import h5py as h5
from . import general_utils as gu__
from . import access as ac__
from . import global_constants as GC

def create_loom_file(folder_path, expr_matrix, barcodes,
                     acc = None, genes = None):
    """Creates a loom file.

    Creates a loom file from a passed expression matrix, list
    of cell barcodes, and list of genes. Also adds the fields
    for the cell-specific universal internal metadata, and
    the HDF5 groups for the author-annotated internal metadata.

    Args:
        folder_path: String. A full path to the folder where
            the new loom file will be stored.
            e.g. "path/to/folder"
            NOT  "path/to/folder/file.loom"
        expr_matrix: The actual expression matrix. Can be a
            numpy.ndarray, or any of these scipy sparse matrices:
                * scipy.sparse.coo_matrix
                * scipy.sparse.csc_matrix
                * scipy.sparse.csr_matrix
            and it must have cells on the columns and
            genes on the rows.
        barcodes: numpy.ndarray. Should be a single dimension
            array holding the cell IDs. Length must equal
            the number of columns in expr_matrix.
        acc: numpy.ndarray. Should be a single dimension
            array holding gene IDs to be stored in Accession.
            These should be unique gene IDs, such as Ensembl
            IDs. Length must equal the number of rows
            in expr_matrix.
        genes: numpy.ndarray. Should be a single dimension
            array holding gene IDs to be stored in Gene.
            These should be human-readable gene IDs and don't
            necessarily have to be unique. Length must equal
            the number of rows in expr_matrix.
    
    Returns: None

    Raises:
        FileNotFoundError: If the passed folder does not exist.
        FileExistsError: If the expr_mat.loom file already exists.
        ValueError: If neither acc nor genes are passed, or if
            acc has non-unique values.
        MemoryError: If loompy.create() fails because of a
            memory overflow error.
    """
    if not os.path.exists(folder_path):
        raise FileNotFoundError(f'{folder_path} doesn\'t exist!')
    file_path = os.path.join(folder_path, 'expr_mat.loom')
    if os.path.exists(file_path):
        raise FileExistsError(f'{file_path} already exists!')
    if acc is None and genes is None:
        raise ValueError('acc and genes cannot both be None!')
    if barcodes.shape != np.unique(barcodes).shape:
        raise ValueError('barcodes cannot have non-unique values!')
    col_attrs = {'CellID': np.array(barcodes)}
    row_attrs = {}
    if acc is not None:
        if acc.shape != np.unique(acc).shape:
            raise ValueError('acc cannot have non-unique values!')
        row_attrs['Accession'] = np.array(acc)
    else:
        row_attrs['Accession'] = np.repeat('-1', expr_matrix.shape[0])
    if genes is not None:
        row_attrs['Gene'] = np.array(genes)
    else:
        row_attrs['Gene'] = np.repeat('-1', expr_matrix.shape[0])
    try:
        lp.create(file_path, expr_matrix, row_attrs, col_attrs)
    except MemoryError as e:
        raise ValueError('Dataset too large to create a loom file!') from e
    with h5.File(file_path, 'r+') as lfile:
        # Add edits to the loom file specification
        lfile.create_group('cell_author_annot')
        lfile.create_group('gene_author_annot')
        lfile['cell_author_annot'].attrs['column_order'] = ''
        lfile['gene_author_annot'].attrs['column_order'] = ''
        cellid_all_missing = (lfile['col_attrs/CellID'][:] == '-1').all()
        lfile['col_attrs/CellID'].attrs['all_missing'] = cellid_all_missing
        if 'Accession' in row_attrs.keys():
            acc_all_missing = (lfile['row_attrs/Accession'][:] == '-1').all()
            lfile['row_attrs/Accession'].attrs['all_missing'] = acc_all_missing
        if 'Gene' in row_attrs.keys():
            gene_all_missing = (lfile['row_attrs/Gene'][:] == '-1').all()
            lfile['row_attrs/Gene'].attrs['all_missing'] = gene_all_missing

def get_expr_matrix_from_cr_triplet(path, prefix, sparse=True):
    """Gets mat, bar, feat from Cell Ranger triplet files.

    Converts the barcodes.tsv, genes/features.tsv, and
    matrix.mtx files outputted by Cell Ranger into:
        1. A scipy sparse CSC matrix holding the expression
           matrix.
        2. A numpy ndarray holding the barcodes
        3. A numpy ndarray holding the features
    However, if there is additional information in the
    barcodes and features files, a pandas dataframe is
    returned instead of a numpy ndarray.

    Expects to receive the path to the folder holding the
    actual barcodes.tsv, features.tsv, matrix.mtx files.
    Can handle gzipped versions (e.g. barcodes.tsv.gz) and
    v2 where the features file is named genes.tsv.

    Args:
        path: String to the folder holding the files
        prefix: String. Prefix appended to the beginning of the
            expected filenames. e.g. Sample1_barcodes.tsv
        sparse: Boolean. If True, a scipy sparse CSC matrix
            will be returned. If False, a numpy ndarray
            will be returned.

    Returns:
        A tuple with 3 members:
            1. A matrix holding the expression
               matrix. Barcodes = Columns; Features = Rows.
               If sparse=True, this will be a scipy CSC
               sparse matrix. Otherwise, will be a numpy
               ndarray.
            2. If the barcodes.tsv file only held the barcodes,
               a numpy ndarray holding the barcodes is returned.
               If there was any additional information, a pandas
               dataframe is returned instead and an alert
               is printed.
            3. If the features.tsv file only held the features,
               a numpy ndarray holding the features is returned.
               If there was any additional information, a pandas
               dataframe is returned instead and an alert
               is printed.
    
    Raises:
        FileNotFoundError: If any of the expected files are not
            found in the folder given by path.
    """
    files = os.listdir(path)
    barcode_file = ''
    feature_file = ''
    matrix_file = ''
    if prefix + 'barcodes.tsv' in files:
        barcode_file = os.path.join(path, prefix + 'barcodes.tsv')
    elif prefix + 'barcodes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix + 'barcodes.tsv.gz'))
        barcode_file = os.path.join(path, prefix + 'barcodes.tsv')
    else:
        raise FileNotFoundError('barcodes file not found! Must be named '
                                '*barcodes.tsv or *barcodes.tsv.gz')

    if prefix + 'features.tsv' in files:
        feature_file = os.path.join(path, prefix + 'features.tsv')
    elif prefix + 'features.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix + 'features.tsv.gz'))
        feature_file = os.path.join(path, prefix + 'features.tsv')
    elif prefix + 'genes.tsv' in files:
        feature_file = os.path.join(path, prefix + 'genes.tsv')
    elif prefix + 'genes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix + 'genes.tsv.gz'))
        feature_file = os.path.join(path, prefix + 'genes.tsv')
    else:
        raise FileNotFoundError('features file not found! Must be named '
                                '*features.tsv, *features.tsv.gz, *genes.tsv, '
                                'or *genes.tsv.gz')

    if prefix + 'matrix.mtx' in files:
        matrix_file = os.path.join(path, prefix + 'matrix.mtx')
    elif prefix + 'matrix.mtx.gz' in files:
        gu__.gunzip(os.path.join(path, prefix + 'matrix.mtx.gz'))
        matrix_file = os.path.join(path, prefix + 'matrix.mtx')
    else:
        raise FileNotFoundError('matrix file not found! Must be named '
                                '*matrix.mtx or *matrix.mtx.gz')

    mat = scipy.io.mmread(matrix_file)
    if not sparse:
        mat = mat.toarray()
    else:
        mat = mat.tocsc()
    features_df = pd.read_csv(feature_file, sep = '\t', header = None)
    if len(features_df.columns) > 1:
        print('ALERT: features has extra columns')
        features = features_df
    else:
        features = np.array(features_df.iloc[:, 0])
    barcodes_df = pd.read_csv(barcode_file, sep = '\t', header = None)
    if len(barcodes_df.columns) > 1:
        print('ALERT: barcodes has extra columns')
        barcodes = barcodes_df
    else:
        barcodes = np.array(barcodes_df.iloc[:, 0])
    return(mat, barcodes, features)

def get_expr_matrix_from_cr_filt_h5(path, sparse=True):
    """Gets mat, bar, feat from Cell Ranger filtered_feature_bc_matrix.

    Converts the filtered_feature_bc_matrix.h5 file outputted
    by Cell Ranger into:
        1. A scipy sparse CSC matrix holding the expression
           matrix.
        2. A numpy ndarray holding the barcodes
        3. A numpy ndarray holding the features
    However, if there is additional information on the
    barcodes and features, a pandas dataframe is
    returned instead of a numpy ndarray.

    Expects to receive the path to the file originally
    named filtered_feature_bc_matrix.h5.

    Args:
        path: String to the file
        sparse: Boolean. If True, a sparse scipy CSC matrix
            will be returned. If False, a numpy ndarray
            will be returned

    Returns:
        A tuple with 3 members:
            1. A matrix holding the expression
               matrix. Barcodes = Columns; Features = Rows.
               If sparse=True, this will be a scipy CSC
               sparse matrix. Otherwise, will be a numpy
               ndarray.
            2. If the barcodes information only held the barcodes,
               a numpy ndarray holding the barcodes is returned.
               If there was any additional information, a pandas
               dataframe is returned instead and an alert
               is printed.
            3. If the features information only held the features,
               a numpy ndarray holding the features is returned.
               If there was any additional information, a pandas
               dataframe is returned instead and an alert
               is printed.
    
    Raises:
        AssertionError: If anything about the HDF5 file structure
            is unexpected.
    """
    with h5.File(path, 'r') as hfile:
        if len(hfile.keys()) != 1:
            raise AssertionError('HDF5 file structure is unusual! '
                                 'You should parse manually.')
        top_group = list(hfile.keys())[0]
        valid_keys = ('barcodes',
                      'data',
                      'indices',
                      'indptr',
                      'shape',
                      'genes',
                      'gene_names',
                      'features',
                      'feature_names')
        for k in hfile[f'{top_group}'].keys():
            if k not in valid_keys:
                raise AssertionError('HDF5 file structure is unusual! '
                                    'You should parse manually.')
        keys = hfile[f'{top_group}'].keys()
        if (('genes' in keys or 'gene_names' in keys)
                and ('features' in keys or 'feature_names' in keys)):
            raise AssertionError('HDF5 file structure is unusual! '
                                 'You should parse manually.')
        elif ('genes' in keys or 'gene_names' in keys):
            feat_name = 'gene'
        elif ('features' in keys or 'feature_names' in keys):
            feat_name = 'feature'
        else:
            print('WARNING: No feature information available!')
            feat_name = None
        mandatory_keys = ('barcodes',
                          'data',
                          'indices',
                          'indptr',
                          'shape')
        for mk in mandatory_keys:
            if mk not in keys:
                raise AssertionError('HDF5 file structure is unusual! '
                                    'You should parse manually.')
        if hfile[f'{top_group}/shape'][1] != hfile[f'{top_group}/barcodes'].shape[0]:
            raise AssertionError('barcodes and shape do not match! '
                                 'You should parse manually.')
        mat = scipy.sparse.csc_matrix((hfile[f'{top_group}/data'][:],
                                       hfile[f'{top_group}/indices'][:],
                                       hfile[f'{top_group}/indptr'][:]),
                                      shape=hfile[f'{top_group}/shape'][:])
        if not sparse:
            mat = mat.toarray()
        bar = hfile[f'{top_group}/barcodes'][:]
        if type(bar[0]) is np.bytes_:
            bar = np.array(list(map(lambda x: x.decode('utf-8'), bar)))
        if feat_name is not None:
            if f'{feat_name}s' in keys and f'{feat_name}_names' in keys:
                col1 = hfile[f'{top_group}/{feat_name}s'][:]
                col2 = hfile[f'{top_group}/{feat_name}_names'][:]
                if col1.shape[0] != col2.shape[0]:
                    raise AssertionError(f'{feat_name}s and {feat_name}_names '
                                          'have different shapes!')
                if type(col1[0]) is np.bytes_:
                    col1 = map(lambda x: x.decode('utf-8'), col1)
                if type(col2[0]) is np.bytes_:
                    col2 = map(lambda x: x.decode('utf-8'), col2)
                # You can pass iterators to pd.DataFrame,
                # so you don't have to convert to an iterable
                feat = pd.DataFrame([col1, col2]).transpose()
            elif f'{feat_name}s' in keys:
                feat = hfile[f'{top_group}/{feat_name}s'][:]
            else:
                feat = hfile[f'{top_group}/{feat_name}_names'][:]
        if feat.shape[0] != hfile[f'{top_group}/shape'][0]:
            raise AssertionError('feature info and shape do not match '
                                 'dimensions! You should parse manually.')
        return(mat, bar, feat)        

def get_expr_matrix_from_csv(path, rows = 'genes', **kwargs):
    """Gets expression matrix, barcodes, features from a csv-like file.

    Converts a csv-like file with a header and rownames into:
        1. A numpy ndarray holding the expression matrix.
        2. A numpy ndarray holding the barcodes
        3. A numpy ndarray holding the features

    Args:
        path: String. Path to the csv-like file.
        rows: String. Must be 'genes' or 'cells'. Indicates the
            rows variable in the file.
        **kwargs: Passed to pandas.read_csv()

    Returns:
        A tuple with 3 members:
            1. A numpy ndarray holding the expression matrix.
               Barcodes = Rows; Features = Columns
            2. A numpy ndarray holding the barcodes
            3. A numpy ndarray holding the features
    
    Raises:
        ValueError: If the rows argument isn't "genes" or "cells"
    """
    df = pd.read_csv(path, **kwargs)
    if rows == 'cells':
        df = df.transpose()
    elif rows != 'genes':
        raise ValueError('rows must be \'genes\' or \'cells\'!')
    mat = np.array(df)
    barcodes = np.array(df.columns)
    features = np.array(df.index)
    return(mat, barcodes, features)

def set_layer(new_mat, name, uuid = None):
    """Sets an expression matrix.

    Sets an expression matrix in the given dataset.

    Args:
        new_mat: numpy 2-D array of ints or floats. The
            expression matrix to be stored in the dataset.
            Must match the dimensions of the matrix under
            lfile['matrix'].
        name: String. The name of the expression matrix to
            write to. If "matrix", will write to the main
            matrix under lfile['matrix']. Will ask for
            confirmation if a layer already exists or if
            the main matrix is being written to.
        uuid: String. The UUID of the desired dataset.

    Returns: None

    Raises:
        ValueError: If new_mat has the wrong shape.
    """
    new_mat = np.array(new_mat)
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if new_mat.shape != lfile['matrix'].shape:
            raise ValueError('new_mat has the wrong shape!')
        if name == 'matrix':
            print('!!! Warning !!!')
            print('You are overwriting the main matrix.')
            if not gu__.get_yes_or_no('Do you want to continue? (y/n): '):
                return
            else:
                lfile['matrix'] = new_mat
        elif name in lfile['layers'].keys():
            print(f'Layer \"{name}\" already exists!')
            if not gu__.get_yes_or_no('Do you want to overwrite? (y/n): '):
                return
            else:
                lfile[f'layers/{name}'] = new_mat
        else:
            lfile.create_dataset(f'layers/{name}',
                                 data = new_mat)

def main():
    pass

if __name__ == '__main__':
    main()