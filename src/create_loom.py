import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import os
import loompy as lp
import numpy as np
import pandas as pd
import scipy
import h5py as h5
import general_utils as gu__
import global_constants as GC

def create_loom_file(file_path, expr_matrix, barcodes, features,
                     feat_acc = True):
    """Creates a loom file.

    Creates a loom file from a passed expression matrix, list
    of cell barcodes, and list of genes. Also adds the fields
    for the cell-specific universal internal metadata, and
    the HDF5 groups for the author-annotated internal metadata.

    Args:
        file_path: String. A full path, including the file
            name, of the new loom file. e.g. path/to/new_file.loom
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
        features: numpy.ndarray. Should be a single dimension
            array holding the gene IDs. Length must equal
            the number of rows in expr_matrix.
        feat_acc: Boolean. True if the features should be
            stored under 'Accession' in the loom file. This
            is appropriate for unique Ensemble IDs for instance.
            If the features are human-readable gene names that
            may not be unique, then they will be stored under
            'Gene', and this argument should be False.
    
    Returns: None
    Raises: None
    """
    if os.path.exists(file_path):
        raise FileExistsError(f'{file_path} already exists!')
    col_attrs = {'CellID': barcodes}
    if feat_acc:
        row_attrs = {'Accession': features}
    else:
        row_attrs = {'Gene': features}
    try:
        lp.create(file_path, expr_matrix, row_attrs, col_attrs)
    except MemoryError:
        print('\nERROR: Dataset too large to create a loom file!\n')
        sys.exit(1)
    with h5.File(file_path, 'r+') as lfile:

        # Add edits to the loom file specification
        # FLAG
        # lfile.create_group("col_attrs/author_annot")
        # lfile.create_group("row_attrs/author_annot")
        lfile.create_group("cell_author_annot")
        lfile.create_group("gene_author_annot")

    # I'm not sure if I want to keep this. This pre-creates the
    # universal internal metadata fields as empty HDF5 datasets.
    # 
    # n_col = lfile['matrix'].shape[1]
    # # Universal Internal Metadata fields in correct order
    # univ_fields = [item[0] for item in sorted(GC._IMU_COLUMN_INDEX.items(),
    #                                           key = lambda x: x[1])]
    # for field in univ_fields:
    #     lfile.create_dataset(f'col_attrs/{field}', shape = (n_col, ))


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
    
    Raises: None
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


def main():
    pass

if __name__ == '__main__':
    main()