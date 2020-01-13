import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import numpy as np
import pandas as pd
import h5py as h5
import access as ac__
import global_constants as GC

# Accidentally wrote these methods twice. Once here, and once
# access.py. The ones in access.py are better, so I'm copying
# them here. Once I vet them and run some tests, I'll actually
# delete these, but for now, I don't want to lose the code.
#
# def get_cell_int_md_univ(uuid):
#     """Gets the cell internal universal metadata for a dataset.
# 
#     For a given dataset, get the cell internal universal
#     metadata as a Pandas DataFrame.
# 
#     Args:
#         uuid: The UUID of the desired dataset.
# 
#     Returns:
#         A Pandas DataFrame containing 
#     """
#     lfile = ac__.get_h5_conn(uuid)
#     new_df = {}
#     for k, v in lfile['col_attrs'].items():
#         if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[1]):
#             new_df[k] = pd.Series(v)
#     idx = new_df['CellID']
#     del new_df['CellID']
#     new_df = pd.DataFrame(new_df)
#     new_df.index = idx
#     lfile.close()
#     return(new_df)
# 
# def get_gene_int_md_univ(uuid):
#     lfile = ac__.get_h5_conn(uuid)
#     new_df = {}
#     for k, v in lfile['row_attrs'].items():
#         if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[0]):
#             new_df[k] = pd.Series(v)
#     if 'Gene' in lfile['row_attrs'].keys():
#         idx = new_df['Gene']
#         del new_df['Gene']
#     elif 'Accession' in lfile['row_attrs'].keys():
#         idx = new_df['Accession']
#         del new_df['Accession']
#     else:
#         raise AssertionError(f'Loom file {uuid} does not have a '
#                               'row_attrs/Gene or a row_attrs/Accession!')
#     new_df = pd.DataFrame(new_df)
#     new_df.index = idx
#     lfile.close()
#     return(new_df)
# 
# def get_cell_int_md_author_annot(uuid):
#     lfile = ac__.get_h5_conn(uuid)
#     new_df = {}
#     for k, v in lfile['col_attrs/author_annot'].items():
#         if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[1]):
#             new_df[k] = pd.Series(v)
#     idx = pd.Series(lfile['col_attrs/CellID'])
#     new_df = pd.DataFrame(new_df)
#     new_df.index = idx
#     lfile.close()
#     return(new_df)
# 
# def get_gene_int_md_author_annot(uuid):
#     lfile = ac__.get_h5_conn(uuid)
#     new_df = {}
#     for k, v in lfile['row_attrs/author_annot'].items():
#         if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[0]):
#             new_df[k] = pd.Series(v)
#     if 'Gene' in lfile['row_attrs'].keys():
#         idx = pd.Series(lfile['row_attrs/Gene'])
#     elif 'Accession' in lfile['row_attrs'].keys():
#         idx = pd.Series(lfile['row_attrs/Accession'])
#     else:
#         raise AssertionError(f'Loom file {uuid} does not have a '
#                               'row_attrs/Gene or a row_attrs/Accession!')
#     new_df = pd.DataFrame(new_df)
#     new_df.index = idx
#     lfile.close()
#     return(new_df)
# 
# def get_cell_int_md(uuid, universal = True):
#     if universal:
#         return(get_cell_int_md_univ(uuid))
#     else:
#         return(get_cell_int_md_author_annot(uuid))
# 
# def get_gene_int_md(uuid, universal = True):
#     if universal:
#         return(get_gene_int_md_univ(uuid))
#     else:
#         return(get_gene_int_md_author_annot(uuid))

def get_cell_int_md_univ(uuid, keep_missing = True):
    """Get the cell universal metadata.

    Get the cell universal metadata from the given dataset
    as a Pandas DataFrame. 

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
    lfile = ac__.get_h5_conn(uuid)
    col_keys = list(lfile['col_attrs'].keys())
    col_data = pd.DataFrame(index = lfile['col_attrs/CellID'])
    for key in col_keys:
        if key == 'author_annot':
            continue
        key_path = 'col_attrs/' + key
        if len(lfile[key_path].shape) == 1:
            if (np.array(lfile[key_path], dtype = '<U1') == '-1').all():
                if keep_missing:
                    col_data[key] = lfile[key_path][:]
            else:
                col_data[key] = lfile[key_path][:]
    if col_data.shape[1] == 0:
        col_data = None
    lfile.close()
    return(col_data)

def get_gene_int_md_univ(uuid, keep_missing = True):
    """Get the gene universal metadata.

    Get the gene universal metadata from the given dataset
    as a Pandas DataFrame. 

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
    lfile = ac__.get_h5_conn(uuid)
    row_keys = list(lfile['row_attrs'].keys())
    if 'Gene' in row_keys:
        row_data = pd.DataFrame(index = lfile['row_attrs/Gene'])
    elif 'Accession' in row_keys:
        row_data = pd.DataFrame(index = lfile['row_attrs/Accession'])
    else:
        row_data = pd.DataFrame()
    for key in row_keys:
        if key == 'author_annot':
            continue
        key_path = 'row_attrs/' + key
        if len(lfile[key_path].shape) == 1:
            if (np.array(lfile[key_path], dtype = '<U1') == '-1').all():
                if keep_missing:
                    row_data[key] = lfile[key_path][:]
            else:
                row_data[key] = lfile[key_path][:]
    if row_data.shape[1] == 0:
        row_data = None
    lfile.close()
    return(row_data)

def get_cell_int_md_author_annot(uuid):
    """Get the cell author-annotated metadata.

    Get the cell author-annotated metadata from the given
    dataset as a Pandas DataFrame.

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as cells
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    lfile = ac__.get_h5_conn(uuid)
    col_keys = list(lfile['col_attrs/author_annot'].keys())
    col_data = pd.DataFrame(index = lfile['col_attrs/CellID'])
    for key in col_keys:
        key_path = 'col_attrs/author_annot/' + key
        if len(lfile[key_path].shape) == 1:
            col_data[key] = lfile[key_path][:]
    if col_data.shape[1] == 0:
        col_data = None
    lfile.close()
    return(col_data)

def get_gene_int_md_author_annot(uuid):
    """Get the gene author-annotated metadata.

    Get the gene author-annotated metadata from the given
    dataset as a Pandas DataFrame.

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as genes
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    lfile = ac__.get_h5_conn(uuid)
    row_keys = list(lfile['row_attrs/author_annot'].keys())
    if 'Gene' in row_keys:
        row_data = pd.DataFrame(index = lfile['row_attrs/Gene'])
    elif 'Accession' in row_keys:
        row_data = pd.DataFrame(index = lfile['row_attrs/Accession'])
    else:
        row_data = pd.DataFrame()
    for key in row_keys:
        key_path = 'row_attrs/author_annot/' + key
        if len(lfile[key_path].shape) == 1:
            row_data[key] = lfile[key_path][:]
    if row_data.shape[1] == 0:
        row_data = None
    lfile.close()
    return(row_data)

def set_cell_int_md_author_annot(uuid, df):
    """Set the cell author-annotated metadata.

    Set the cell author-annotated metadata for a given dataset
    from a Pandas DataFrame. Confirmation must be given
    before overwriting existing columns. Additionally, the method
    will attempt to undo all of its edits and exit gracefully
    if any error occurs. Use this method carefully to avoid
    losing existing metadata.

    Args:
        uuid: String. UUID of the desired dataset.
        df: Pandas DataFrame. Each column should be a field
            from the cell author-annotated metadata.
    
    Returns: None

    Raises:
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[1]:
        raise AssertionError('df has the wrong number of rows!')
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['col_attrs/author_annot'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'col_attrs/author_annot/{col}'])
                    del lfile[f'col_attrs/author_annot/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'col_attrs/author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['col_attrs/author_annot'].keys():
                        del lfile[f'col_attrs/author_annot/{col_del}']
                    dset = lfile.create_dataset(f'col_attrs/author_annot/{col_del}',
                                                (lfile['matrix'].shape[1], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))
    lfile.close()

def set_gene_int_md_author_annot(uuid, df):
    """Set the gene author-annotated metadata.

    Set the gene author-annotated metadata for a given dataset
    from a Pandas DataFrame. Confirmation must be given
    before overwriting existing columns. Additionally, the method
    will attempt to undo all of its edits and exit gracefully
    if any error occurs. Use this method carefully to avoid
    losing existing metadata.

    Args:
        uuid: String. UUID of the desired dataset.
        df: Pandas DataFrame. Each column should be a field
            from the cell author-annotated metadata.
    
    Returns: None

    Raises:
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[0]:
        raise AssertionError('df has the wrong number of rows!')
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['row_attrs/author_annot'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'row_attrs/author_annot/{col}'])
                    del lfile[f'col_attrs/author_annot/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'row_attrs/author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['row_attrs/author_annot'].keys():
                        del lfile[f'row_attrs/author_annot/{col_del}']
                    dset = lfile.create_dataset(f'row_attrs/author_annot/{col_del}',
                                                (lfile['matrix'].shape[1], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))
    lfile.close()

def set_cell_int_md_univ(uuid, df):
    """Set the cell universal metadata.

    Set the cell universal metadata for a given dataset
    from a Pandas DataFrame. Confirmation must be given
    before overwriting existing columns. Additionally, the method
    will attempt to undo all of its edits and exit gracefully
    if any error occurs. Use this method carefully to avoid
    losing existing metadata.

    Args:
        uuid: String. UUID of the desired dataset.
        df: Pandas DataFrame. Each column should be an appropriately
            named and typed field from the universal metadata
            schema. There is no validation against the universal
            metadata schema.
    
    Returns: None

    Raises:
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[1]:
        raise ValueError('df has the wrong number of rows!')
    test_cols = np.array(df.columns)
    if not np.isin(test_cols, np.array(GC._IMU_COL_COLUMN_INDEX.keys())).all():
        raise ValueError('df has invalid columns!')
    if np.unique(test_cols).shape[0] != df.columns.shape[0]:
        raise ValueError('df has non-unique columns!')
    del test_cols
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['col_attrs'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'col_attrs/{col}'])
                    del lfile[f'col_attrs/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'col_attrs/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'col_attrs/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'col_attrs/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['col_attrs/'].keys():
                        del lfile[f'col_attrs/{col_del}']
                    dset = lfile.create_dataset(f'col_attrs/{col_del}',
                                                (lfile['matrix'].shape[1], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))
    lfile.close()

def set_gene_int_md_univ(uuid, df):
    """Set the gene universal metadata.

    Set the gene universal metadata for a given dataset
    from a Pandas DataFrame. Confirmation must be given
    before overwriting existing columns. Additionally, the method
    will attempt to undo all of its edits and exit gracefully
    if any error occurs. Use this method carefully to avoid
    losing existing metadata.

    Args:
        uuid: String. UUID of the desired dataset.
        df: Pandas DataFrame. Each column should be an appropriately
            named and typed field from the universal metadata
            schema. There is no validation against the universal
            metadata schema.
    
    Returns: None

    Raises:
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[0]:
        raise ValueError('df has the wrong number of rows!')
    test_cols = np.array(df.columns)
    if not np.isin(test_cols, np.array(GC._IMU_ROW_COLUMN_INDEX.keys())).all():
        raise ValueError('df has invalid columns!')
    if np.unique(test_cols).shape[0] != df.columns.shape[0]:
        raise ValueError('df has non-unique columns!')
    del test_cols
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['row_attrs'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'row_attrs/{col}'])
                    del lfile[f'row_attrs/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'row_attrs/{col}',
                                            (lfile['matrix'].shape[0], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'row_attrs/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'row_attrs/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['row_attrs/'].keys():
                        del lfile[f'row_attrs/{col_del}']
                    dset = lfile.create_dataset(f'row_attrs/{col_del}',
                                                (lfile['matrix'].shape[0], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))
    lfile.close()


def main():
    pass

if __name__ == '__main__':
    main()