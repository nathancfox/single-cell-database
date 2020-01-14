import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
# Used for section below
# import re
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

def construct_batch(df, columns):
    """Construct a universal internal metadata batch field.

    Given a dataframe and a list of columns in that dataframe,
    concatenate the values to construct a 'batch' column
    for the cell internal metadata.

    Args:
        df: Pandas DataFrame. Contains columns that will be
            parts of the batch column.
        columns: List of strings. List of column names,
            in order, that will be concatenated to
            form the 'batch' column.

    Returns:
        A tuple with 2 members:
          1. A Pandas Series object containing the new 'batch'
             column for the cell internal universal metadata.
          2. A single string with the column names concatenated
             to describe the new batch column values.

    Raises: None
    """
    df = df[columns]
    out = []
    for i in range(df.shape[0]):
        out.append('|'.join(df.iloc[i].apply(str)))
    out = pd.Series(out)
    return((out, '|'.join(columns)))

# These methods were built to support a MAPPING schema
# originally described in blog post "Nathan: Jan 6 - Jan 10"
# I decided not to use the MAPPING schema, at least at
# first and so all of these are UNTESTED!!! and not
# being used right now.
# 
# def construct_default_mapping(column, zero_value = None,
#                               batch = False):
#     """Construct a mapping for a universal metadata field.

#     Given a list of values, convert the list to a list
#     of integers, mapped to the original set of unique
#     values. Also return the array of strings to be
#     set as the HDF5 attribute for this HDF5 dataset.
#     See the internal universal metadata specification
#     for details.

#     Args:
#         column: List-like object. List of values to be
#             converted to a mapping. Must be convertable
#             to a Pandas Series.
#         zero_value: Any of the values of column. If not
#             None, this value will be mapped to 0 instead
#             of a regular mapping integer.
#         batch: If True, handles the case where column is
#             the 'batch' column and there are multiple
#             '|'-delimited fields inside the mapping.
    
#     Returns:
#         A 2-member tuple. The first member is the original list
#         with the values replaced by their corresponding integers.
#         The second member is an array of mappings from the integers
#         to their descriptions.

#     Raises: None
#     """
#     column = pd.Series(column) 
#     if batch:
#         raise AssertionError('batch is not built yet!')
#     # TODO(Nathan Fox): Convert NaN etc to "-1"
#     mapping_dict = {}
#     counter = 1
#     for value in column.unique():
#         if zero_value is not None:
#             if value == zero_value:
#                 mapping_dict[str(value)] = 0
#                 continue
#         mapping_dict[str(value)] = counter
#         counter += 1
#     mapping = column.map(mapping_dict)
#     rev_mapping_dict = {}
#     for k, v in mapping_dict:
#         rev_mapping_dict[str(v)] = k
#     mapping_key = []
#     for k in sorted(list(rev_mapping_dict.keys())):
#         mapping_key.append(f'{k}:{rev_mapping_dict[k]}')
#     mapping_key = np.array(mapping_key)
#     return((mapping, mapping_key)) 

# def construct_interactive_mapping(column, batch = False):
#     """Run an interactive constructor for a universal metadata mapping.

#     Run an interactive input application for constructing an
#     internal universal metadata mapping.

#     Args:
#         column: List-like object. List of values to be
#             converted to a mapping. Must be convertable
#             to a Pandas Series.
#         batch: If True, handles the case where column is
#             the 'batch' column and there are multiple
#             '|'-delimited fields inside the mapping.
    
#     Returns:
#         A 2-member tuple. The first member is the original list
#         with the values replaced by their corresponding integers.
#         The second member is an array of mappings from the integers
#         to their descriptions.

#     Raises: None
#     """
#     column = pd.Series(column) 
#     if batch:
#         raise AssertionError('batch is not built yet!')
#     # TODO(Nathan Fox): Convert NaN etc to "-1"
#     mapping_dict = {}
#     counter = 1
#     for value in column.unique():
#         if zero_value is not None:
#             if value == zero_value:
#                 mapping_dict[str(value)] = 0
#                 continue
#         mapping_dict[str(value)] = counter
#         counter += 1
#     mapping = column.map(mapping_dict)
#     rev_mapping_dict = {}
#     for k, v in mapping_dict:
#         rev_mapping_dict[str(v)] = k
#     intro_text = (
#                     '\n'
#                     'Construct a Mapping\n'
#                     '-------------------\n'
#                     '\n'
#                     'For each unique value in the Series, map it to\n'
#                     'an integer and store its description in a mapping.\n'
#                     'Each mapping value will be stored as a string of\n'
#                     'the format \"INTEGER:DESCRIPTION\", where the mapped\n' 
#                     'integer is followed by a colon, then a string of any\n'
#                     'length describing the value. Each prompt will ask you\n'
#                     'for the DESCRIPTION text. The text may not have\n'
#                     'a vertical bar character \"|\" in it. If the substring\n'
#                     '\"{name}\" is in the DESCRIPTION, it will be replaced\n'
#                     'with the name of the value.\n'
#                     '\n'
#                     'To execute one of the following options, enter\n'
#                     'the command in any of the prompts.\n'
#                     '    QUIT    : Abort this mapping construction\n'
#                     '    RESTART : Start this mapping construction over\n'
#                     '    HELP    : Get instructions for input\n'
#                     '    LIST    : Get a list of already entered mappings\n'
#     ) 
#     print(intro_text)
#     input('Press Enter to continue...')
#     print()
#     loop = True
#     unique_vals = list(mapping_dict.keys())
#     val_idx = 0
#     # Not sure if it's better for mapping_key to be a dict
#     # or a list.
#     mapping_key = ['FILL_ME_IN' for i in range(len(unique_vals))]
#     # mapping_key = {val: 'FILL_ME_IN' for val in unique_vals}
#     print_uniq_vals = input(f'There are {len(unique_vals)}.\n'
#                              'Do you want to view them before '
#                              'starting? (y/n) : ')[0].lower()
#     if print_uniq_vals == 'y':
#         print()
#         print('Unique Values')
#         print('-------------')
#         for val in unique_vals:
#             print(f'  {val}')
#         print()
#     while loop:
#         val_name = unique_vals[val_idx]
#         print(f'Value: {val_name}')
#         print('Description:')
#         user_input = input('> ')
#         if user_input == 'QUIT':
#             print('\nExiting...')
#             return
#         elif user_input == 'RESTART':
#             print('\nStarting over...')
#             print(intro_text)
#             loop = True
#             unique_vals = list(mapping_dict.keys())
#             val_idx = 0
#             mapping_key = ['FILL_ME_IN' for i in range(len(unique_vals))]
#             print_uniq_vals = input(f'There are {len(unique_vals)}.\n'
#                                     'Do you want to view them before '
#                                     'starting? (y/n) : ')[0].lower()
#             if print_uniq_vals == 'y':
#                 print()
#                 print('Unique Values')
#                 print('-------------')
#                 for val in unique_vals:
#                     print(f'  {val}')
#                 print()
#             continue
#         elif user_input == 'HELP':
#             help_text = (
#                             'Description Format:\n'
#                             '    1. No \"|\" characters allowed.\n'
#                             '    2. \"{name}\" will be replaced with the value name.\n'
#                             '\n'
#                             'Commands:\n'
#                             '    QUIT    : Abort this mapping construction\n'
#                             '    RESTART : Start this mapping construction over\n'
#                             '    HELP    : Get this help page\n'
#                             '    LIST    : Get a list of already entered mappings\n'
#             )
#             print()
#             print(help_text)
#             print()
#             continue
#         elif user_input == 'LIST':
#             print()
#             print('Values Already Entered')
#             print('----------------------')
#             max_width = max(list(map(lambda x: 0 if x == 'FILL_ME_IN' else len(x),
#                                      mapping_key)))
#             for val in mapping_key:
#                 if val == 'FILL_ME_IN':
#                     continue
#                 val = val.split(':')
#                 print(f'  {rev_mapping_dict[val[0]]:{max_width}s} : {''.join(val[1:])}')
#             print()
#             continue
#         else:
#             if user_input == '':
#                 if input('Are you sure you want to leave this '
#                          'DESCRIPTION empty? (y/n) : ')[0].lower() == 'n':
#                     print()
#                     continue
#             if '|' in user_input:
#                 print()
#                 print('ERROR: You cannot have a \"|\" character in your DESCRIPTION!')
#                 print()
#                 continue
#             new_desc = re.sub(r'{name}', val_name, user_input)
#             new_desc = f'{mapping_dict[val_name]}:{new_desc}'
#             mapping_key[val_idx] = new_desc
#             val_idx += 1
#             if val_idx >= len(unique_vals):
#                 loop = False
#             print()
#     if 'FILL_ME_IN' in mapping_key:
#         raise AssertionError('At least one value did not get a DESCRIPTION!')
#     # if 'FILL_ME_IN' in mapping_key.values():
#     #     raise AssertionError('At least one value did not get a DESCRIPTION!')
#     return((mapping, mapping_key))
                
# def construct_mapping(column, zero_value = None, mode = 'default')
#     """Construct a mapping for an internal universal metadata column.

#     Acts as a router method to one of the two helper methods
#     with appropriate arguments selected.
#     """
#     if mode == 'default':
#         return(construct_default_mapping(column, zero_value = zero_value))
#     elif mode == 'interactive':
#         return(construct_interactive_mapping(column))
#     elif mode == 'batch_default':
#         return(construct_default_mapping(column, zero_value = zero_value,
#                                          batch = True))
#     elif mode == 'batch_interactive':
#         return(construct_interactive_mapping(column, batch = True))
#     else:
#         raise ValueError('mode can only be one of [\"default\", '
#                          '\"interactive\", \"batch_default\", '
#                          '\"batch_interactive\"]')

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
        If there is a "batch_key" attribute for the batch column
        and the batch column isn't a missing column, then a
        2-member tuple is returned:
          1. Pandas DataFrame with the same number of rows as cells
             in the dataset and where each column is a universal
             internal metadata field.
          2. String. The batch_key.
        If not, just the Pandas DataFrame is returned, not in a tuple.

    Raises: None
    """
    lfile = ac__.get_h5_conn(uuid)
    col_keys = list(lfile['col_attrs'].keys())
    col_data = pd.DataFrame(index = np.array(lfile['col_attrs/CellID']))
    for key in col_keys:
        # FLAG
        # if key == 'author_annot' or key == 'CellID':
        if key == 'CellID':
            continue
        key_path = 'col_attrs/' + key
        if len(lfile[key_path].shape) == 1:
            if (np.array(list(map(str, list(lfile[key_path])))) == '-1').all():
                if keep_missing:
                    col_data[key] = lfile[key_path][:]
            else:
                col_data[key] = lfile[key_path][:]
    if col_data.shape[1] == 0:
        col_data = None
        return(col_data)
    else:
        column_order = sorted(col_data.columns,
                            key = lambda x: GC._IMU_CELL_COLUMN_INDEX[x])
        col_data = col_data[column_order]
        batch_key = None
        if ('batch' in col_data.columns
            and 'batch_key' in lfile['col_attrs/batch'].attrs.keys()):
            batch_key = lfile['col_attrs/batch'].attrs['batch_key']
        lfile.close()
        if batch_key is not None:
            return((col_data, batch_key))
        else:
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
        row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Gene']))
        skip_col = 'Gene'
    elif 'Accession' in row_keys:
        row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Accession']))
        skip_col = 'Accession'
    else:
        row_data = pd.DataFrame()
    for key in row_keys:
        # FLAG
        # if key == 'author_annot' or key == skip_col:
        if key == skip_col:
            continue
        key_path = 'row_attrs/' + key
        if len(lfile[key_path].shape) == 1:
            if (np.array(list(map(str, list(lfile[key_path])))) == '-1').all():
                if keep_missing:
                    row_data[key] = lfile[key_path][:]
            else:
                row_data[key] = lfile[key_path][:]
    if row_data.shape[1] == 0:
        row_data = None
    else:
        column_order = sorted(row_data.columns,
                              key = lambda x: GC._IMU_GENE_COLUMN_INDEX[x])
        if skip_col == 'Gene' and 'Accession' in row_keys:
            column_order = ['Accession'] + column_order
        elif skip_col == 'Accession' and 'Gene' in row_keys:
            column_order = ['Gene'] + column_order
        row_data = row_data[column_order]
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
    # FLAG
    # col_keys = list(lfile['col_attrs/author_annot'].keys())
    col_keys = list(lfile['cell_author_annot'].keys())
    col_data = pd.DataFrame(index = np.array(lfile['col_attrs/CellID']))
    for key in col_keys:
        # FLAG
        # key_path = 'col_attrs/author_annot/' + key
        key_path = 'cell_author_annot/' + key
        if len(lfile[key_path].shape) == 1:
            col_data[key] = lfile[key_path][:]
    if col_data.shape[1] == 0:
        col_data = None
    else:
        # FLAG
        # col_order = lfile['col_attrs/author_annot'].attrs['column_order'].split('|')
        col_order = lfile['cell_author_annot'].attrs['column_order'].split('|')
        col_data = col_data[col_order]
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
    # FLAG
    # row_keys = list(lfile['row_attrs/author_annot'].keys())
    row_keys = list(lfile['gene_author_annot'].keys())
    if 'Gene' in row_keys:
        row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Gene']))
    elif 'Accession' in row_keys:
        row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Accession']))
    else:
        row_data = pd.DataFrame()
    for key in row_keys:
        # FLAG
        # key_path = 'row_attrs/author_annot/' + key
        key_path = 'gene_author_annot/' + key
        if len(lfile[key_path].shape) == 1:
            row_data[key] = lfile[key_path][:]
    if row_data.shape[1] == 0:
        row_data = None
    else:
        # FLAG
        # col_order = lfile['row_attrs/author_annot'].attrs['column_order'].split('|')
        col_order = lfile['gene_author_annot'].attrs['column_order'].split('|')
        row_data = row_data[col_order]
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
            # FLAG
            # if col in lfile['col_attrs/author_annot'].keys():
            if col in lfile['cell_author_annot'].keys():
                overwrite_col = input(f'Column \"{col}\" already exists!\n'
                                       'Overwrite? (y/n): ')[0].lower()
                if overwrite_col == 'n':
                    continue
                else:
                    # FLAG
                    # columns_written[col] = pd.Series(lfile[f'col_attrs/author_annot/{col}'])
                    columns_written[col] = pd.Series(lfile[f'cell_author_annot/{col}'])
                    # FLAG
                    # del lfile[f'col_attrs/author_annot/{col}']
                    del lfile[f'cell_author_annot/{col}']
            if df[col].dtype == object:
                # FLAG
                # dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                #                             (lfile['matrix'].shape[1], ),
                #                             dtype = h5.string_dtype())
                dset = lfile.create_dataset(f'cell_author_annot/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                if col not in columns_written.keys():
                    columns_written[col] = None
                dset[:] = df[col]
            else:
                # FLAG
                # dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                #                             data = df[col])
                dset = lfile.create_dataset(f'cell_author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    # FLAG
                    # del lfile[f'col_attrs/author_annot/{col_del}']
                    del lfile[f'cell_author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    # FLAG
                    # if col_del in lfile['col_attrs/author_annot'].keys():
                    #     del lfile[f'col_attrs/author_annot/{col_del}']
                    # dset = lfile.create_dataset(f'col_attrs/author_annot/{col_del}',
                    #                             (lfile['matrix'].shape[1], ),
                    #                             dtype = h5.string_dtype())
                    if col_del in lfile['cell_author_annot'].keys():
                        del lfile[f'cell_author_annot/{col_del}']
                    dset = lfile.create_dataset(f'cell_author_annot/{col_del}',
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
    # FLAG
    # lfile['col_attrs/author_annot'].attrs['column_order'] = '|'.join(df.columns)
    lfile['cell_author_annot'].attrs['column_order'] = '|'.join(df.columns)
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
            # FLAG
            # if col in lfile['row_attrs/author_annot'].keys():
            if col in lfile['gene_author_annot'].keys():
                overwrite_col = input(f'Column \"{col}\" already exists!\n'
                                       'Overwrite? (y/n): ')[0].lower()
                if overwrite_col == 'n':
                    continue
                else:
                    # FLAG
                    # columns_written[col] = pd.Series(lfile[f'row_attrs/author_annot/{col}'])
                    columns_written[col] = pd.Series(lfile[f'gene_author_annot/{col}'])
                    # FLAG
                    # del lfile[f'col_attrs/author_annot/{col}']
                    del lfile[f'cell_author_annot/{col}']
            if df[col].dtype == object:
                # FLAG
                # dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                #                             (lfile['matrix'].shape[0], ),
                #                             dtype = h5.string_dtype())
                dset = lfile.create_dataset(f'gene_author_annot/{col}',
                                            (lfile['matrix'].shape[0], ),
                                            dtype = h5.string_dtype())
                if col not in columns_written.keys():
                    columns_written[col] = None
                dset[:] = df[col]
            else:
                # FLAG
                # dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                #                             data = df[col])
                dset = lfile.create_dataset(f'gene_author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    # FLAG
                    # del lfile[f'row_attrs/author_annot/{col_del}']
                    del lfile[f'gene_author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    # FLAG
                    # if col_del in lfile['row_attrs/author_annot'].keys():
                    #     del lfile[f'row_attrs/author_annot/{col_del}']
                    # dset = lfile.create_dataset(f'row_attrs/author_annot/{col_del}',
                    #                             (lfile['matrix'].shape[1], ),
                    #                             dtype = h5.string_dtype())
                    if col_del in lfile['gene_author_annot'].keys():
                        del lfile[f'gene_author_annot/{col_del}']
                    dset = lfile.create_dataset(f'gene_author_annot/{col_del}',
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
    # FLAG
    # lfile['row_attrs/author_annot'].attrs['column_order'] = '|'.join(df.columns)
    lfile['gene_author_annot'].attrs['column_order'] = '|'.join(df.columns)
    lfile.close()

def set_cell_int_md_univ(uuid, df, batch_key):
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
        batch_key: String. The concatenated column names to be
            stored as an HDF5 attribute on the batch universal
            internal metadata field.
            e.g. "Well|SeqRun|Plate"
    
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
    if not np.isin(test_cols, np.array(list(GC._IMU_CELL_COLUMN_INDEX.keys()))).all():
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
                overwrite_col = input(f'Column \"{col}\" already exists!\n'
                                       'Overwrite? (y/n): ')[0].lower()
                if overwrite_col == 'n':
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
            if col == 'batch':
                dset.attrs['batch_key'] = batch_key
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
    if not np.isin(test_cols, np.array(list(GC._IMU_GENE_COLUMN_INDEX.keys()))).all():
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
                overwrite_col = input(f'Column \"{col}\" already exists!\n'
                                       'Overwrite? (y/n): ')[0].lower()
                if overwrite_col == 'n':
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