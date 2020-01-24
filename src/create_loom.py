"""Functions to parse expression data into a new loom file.

The main function is create_loom_file(). This creates a new
loom file for a new entry from a standard set of arguments.
All other functions exist to create that standard set of
arguments from a variety of expression data formats.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
import sys
sys.path.append('/home/scdb_codebase/single_cell_database/src')
import os
import loompy as lp
import numpy as np
import pandas as pd
import scipy
import h5py as h5
import general_utils as gu__
import access as ac__
import global_constants as GC

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
        feat_acc: Boolean. True if the features should be
            stored under 'Accession' in the loom file. This
            is appropriate for unique Ensembl IDs for instance.
            If the features are human-readable gene names that
            may not be unique, then they will be stored under
            'Gene', and this argument should be False.
    
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
    col_attrs = {'CellID': barcodes}
    row_attrs = {}
    if acc is not None:
        if acc.shape != np.unique(acc).shape:
            raise ValueError('acc cannot have non-unique values!')
        row_attrs['Accession'] = acc
    if genes is not None:
        row_attrs['Gene'] = genes
    try:
        lp.create(file_path, expr_matrix, row_attrs, col_attrs)
    except MemoryError:
        print('\nERROR: Dataset too large to create a loom file!\n')
        sys.exit(1)
    with h5.File(file_path, 'r+') as lfile:
        # Add edits to the loom file specification
        lfile.create_group('cell_author_annot')
        lfile.create_group('gene_author_annot')

def get_expr_matrix_from_cellranger(path, prefix):
    """Gets expression matrix, barcodes, features from Cell Ranger.

    Converts the matrices outputted by Cell Ranger into:
        1. A scipy sparse COO matrix holding the expression
           matrix.
        2. A numpy ndarray holding the barcodes
        3. A numpy ndarray holding the features
    Expects to receive the path to the folder holding the
    actual barcodes.tsv, features.tsv, matrix.mtx files.
    Can handle gzipped versions (e.g. barcodes.tsv.gz) and
    v2 where the features file is named genes.tsv.

    Args:
        path: String to the folder holding the files
        prefix: String. Prefix appended to the beginning of the
            expected filenames. e.g. Sample1_barcodes.tsv

    Returns:
        A tuple with 3 members:
            1. A scipy sparse COO matrix holding the expression
               matrix. Barcodes = Rows; Features = Columns
            2. A numpy ndarray holding the barcodes
            3. A numpy ndarray holding the features
    
    Raises:
        FileNotFoundError: If any of the expected files are not
            found in the folder given by path.
    """
    files = os.listdir(path)
    barcode_file = ''
    feature_file = ''
    matrix_file = ''
    if prefix + 'barcodes.tsv' in files:
        barcode_file = os.path.join(path, prefix, 'barcodes.tsv')
    elif prefix + 'barcodes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix, 'barcodes.tsv.gz'))
        barcode_file = os.path.join(path, prefix, 'barcodes.tsv')
    else:
        raise FileNotFoundError('barcodes file not found! Must be named '
                                '*barcodes.tsv or *barcodes.tsv.gz')

    if prefix + 'features.tsv' in files:
        feature_file = os.path.join(path, prefix, 'features.tsv')
    elif prefix + 'features.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix, 'features.tsv.gz'))
        feature_file = os.path.join(path, prefix, 'features.tsv')
    elif prefix + 'genes.tsv' in files:
        feature_file = os.path.join(path, prefix, 'genes.tsv')
    elif prefix + 'genes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, prefix, 'genes.tsv.gz'))
        feature_file = os.path.join(path, prefix, 'genes.tsv')
    else:
        raise FileNotFoundError('features file not found! Must be named '
                                '*features.tsv, *features.tsv.gz, *genes.tsv, '
                                'or *genes.tsv.gz')

    if prefix + 'matrix.mtx' in files:
        matrix_file = os.path.join(path, prefix, 'matrix.mtx')
    elif prefix + 'matrix.mtx.gz' in files:
        gu__.gunzip(os.path.join(path, prefix, 'matrix.mtx.gz'))
        matrix_file = os.path.join(path, prefix, 'matrix.mtx')
    else:
        raise FileNotFoundError('matrix file not found! Must be named '
                                '*matrix.mtx or *matrix.mtx.gz')

    mat = scipy.io.mmread(matrix_file)
    features_df = pd.read_csv(feature_file, sep = '\t', header = None)
    features = np.array(features_df.iloc[:, 0])
    barcodes_df = pd.read_csv(barcode_file, sep = '\t', header = None)
    barcodes = np.array(barcodes_df.iloc[:, 0])
    return(mat, barcodes, features)

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
    Can be accessed via a UUID lookup from the external
    metadata, or by a direct filepath. Cannot pass
    both uuid and filepath.

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
            dset = lfile.create_dataset(f'layers/{name}',
                                        data = new_mat)

def main():
    pass

if __name__ == '__main__':
    main()