import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import os
import loompy as lp
import numpy as np
import pandas as pd
import scipy
import utils.general_utils as gu__

def create_loom_file(file_path, expr_matrix, barcodes, features):
    """Creates a loom file.

    Args:
        file_path: String. A full path, including the file
            name, of the new loom file. e.g. path/to/new_file.loom
        expr_matrix: The actual expression matrix. Can be a
            numpy.ndarray, or any of these scipy sparse matrices:
                * scipy.sparse.coo_matrix
                * scipy.sparse.csc_matrix
                * scipy.sparse.csr_matrix
            and it must have cells on the rows and genes
            on the columns.
        barcodes: numpy.ndarray. Should be a single dimension
            array holding the cell IDs. Length must equal
            the number of rows in expr_matrix.
        features: numpy.ndarray. Should be a single dimension
            array holding the gene IDs. Length must equal
            the number of columns in expr_matrix.
    
    Returns: None
    Raises: None
    """
    if os.path.exists(file_path):
        raise FileExistsError(f'{file_path} already exists!')
    row_attrs = {'barcode': barcodes}
    col_attrs = {'feature': features}
    try:
        lp.create(file_path, expr_matrix, row_attrs, col_attrs)
    except MemoryError:
        print('\nERROR: Dataset too large to create a loom file!\n')
        sys.exit(1)


def get_expr_matrix_from_cellranger(path):
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
    if 'barcodes.tsv' in files:
        barcode_file = os.path.join(path, 'barcodes.tsv')
    elif 'barcodes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, 'barcodes.tsv.gz'))
        barcode_file = os.path.join(path, 'barcodes.tsv')
    else:
        raise FileNotFoundError('barcodes file not found! Must be named '
                                'barcodes.tsv or barcodes.tsv.gz')

    if 'features.tsv' in files:
        feature_file = os.path.join(path, 'features.tsv')
    elif 'features.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, 'features.tsv.gz'))
        feature_file = os.path.join(path, 'features.tsv')
    elif 'genes.tsv' in files:
        feature_file = os.path.join(path, 'genes.tsv')
    elif 'genes.tsv.gz' in files:
        gu__.gunzip(os.path.join(path, 'genes.tsv.gz'))
        feature_file = os.path.join(path, 'genes.tsv')
    else:
        raise FileNotFoundError('features file not found! Must be named '
                                'features.tsv, features.tsv.gz, genes.tsv, '
                                'or genes.tsv.gz')

    if 'matrix.mtx' in files:
        matrix_file = os.path.join(path, 'matrix.mtx')
    elif 'matrix.mtx.gz' in files:
        gu__.gunzip(os.path.join(path, 'matrix.mtx.gz'))
        matrix_file = os.path.join(path, 'matrix.mtx')
    else:
        raise FileNotFoundError('matrix file not found! Must be named '
                                'matrix.mtx or matrix.mtx.gz')

    mat = scipy.io.mmread(matrix_file).transpose()
    features_df = pd.read_csv(feature_file, sep = '\t', header = None)
    features = np.array(features_df.iloc[:, 0])
    barcodes_df = pd.read_csv(barcode_file, sep = '\t', header = None)
    barcodes = np.array(barcodes_df.iloc[:, 0])

    return(mat, barcodes, features)



def main():
    pass

if __name__ == '__main__':
    main()