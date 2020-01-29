"""Functions to manage the internal metadata.

Functions to read, write, and manage the internal metadata
of the database. This includes both cell-specific and
gene-specific universal and author-annotated internal metadata.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import numpy as np
import pandas as pd
import h5py as h5
import re
from . import access as ac__
from . import global_constants as GC
from . import general_utils as gu__

# These methods were built to support a MAPPING schema
# originally described in blog post "Nathan: Jan 6 - Jan 10"
# I decided not to use the MAPPING schema, at least at
# first and so all of these are UNTESTED!!! and not
# being used right now.
# 
# def construct_default_mapping(column, zero_value = None,
#                               batch = False):
#     """Construct a mapping for a universal metadata field.
# 
#     Given a list of values, convert the list to a list
#     of integers, mapped to the original set of unique
#     values. Also return the array of strings to be
#     set as the HDF5 attribute for this HDF5 dataset.
#     See the internal universal metadata specification
#     for details.
# 
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
#
#     Returns:
#         A 2-member tuple. The first member is the original list
#         with the values replaced by their corresponding integers.
#         The second member is an array of mappings from the integers
#         to their descriptions.
#
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
#
# def construct_interactive_mapping(column, batch = False):
#     """Run an interactive constructor for a universal metadata mapping.
#
#     Run an interactive input application for constructing an
#     internal universal metadata mapping.
#
#     Args:
#         column: List-like object. List of values to be
#             converted to a mapping. Must be convertable
#             to a Pandas Series.
#         batch: If True, handles the case where column is
#             the 'batch' column and there are multiple
#             '|'-delimited fields inside the mapping.
#    
#     Returns:
#         A 2-member tuple. The first member is the original list
#         with the values replaced by their corresponding integers.
#         The second member is an array of mappings from the integers
#         to their descriptions.
#
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
#                
# def construct_mapping(column, zero_value = None, mode = 'default')
#     """Construct a mapping for an internal universal metadata column.
#
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

def construct_batch(df, columns):
    """Constructs a universal internal metadata batch field.

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
    out.index = df.index
    for i, col in enumerate(columns):
        if '/' in col:
            columns[i] = re.sub(r'/', r'_', col)
            print('WARNING: internal_metadata.construct_batch()\n'
                  '  \"/\" characters are not allowed in column names!\n'
                 f'  \"{col}\" was changed to \"{columns[i]}\".\n')
    return((out, '|'.join(columns)))

def get_cell_batch_key(uuid):
    """Get the batch_key for the cell universal metadata.

    The cell universal metadata schema says that the "batch"
    column has to have a string stored with it in an HDF5
    attribute called "batch_key". This method gets that
    batch_key string.

    Args:
        uuid: String. The UUID of the desired dataset.

    Returns:
        The "batch_key" string stored in the dataset.
    
    Raises:
        AssertionError: If the "batch" column or the
            "batch_key" attribute doesn't exist.
    """
    with ac__.get_h5_conn(uuid) as lfile:
        if 'batch' not in lfile['col_attrs/'].keys():
            raise AssertionError('Cell universal metadata doesn\'t have '
                                 'a \"batch\" column!')
        if 'batch_key' not in lfile['col_attrs/batch'].attrs.keys():
            raise AssertionError('Cell universal metadata \"batch\" '
                                 'column doesn\'t have a \"batch_key\"'
                                 'attribute!')
        batch_key = lfile['col_attrs/batch'].attrs['batch_key']
        return(batch_key)

def get_cell_int_md_univ(uuid, keep_missing = True):
    """Gets the cell universal metadata.

    Get the cell universal metadata from the given dataset
    as a Pandas DataFrame. 

    Args:
        uuid: String. UUID of the desired dataset
        keep_missing: Boolean. If True, missing columns will
            be retained as columns of all -1 or "-1".
            If False, they are dropped from the returned DataFrame.
    
    Returns:
        Pandas DataFrame with the same number of rows as cells
        in the dataset and where each column is a universal
        internal metadata field.

    Raises: None
    """
    with ac__.get_h5_conn(uuid) as lfile:
        col_keys = list(lfile['col_attrs'].keys())
        col_data = pd.DataFrame(index = np.array(lfile['col_attrs/CellID']))
        for key in col_keys:
            if key == 'CellID':
                continue
            key_path = 'col_attrs/' + key
            if len(lfile[key_path].shape) == 1:
                if (np.array(list(map(str, list(lfile[key_path])))) == '-1').all():
                    if keep_missing:
                        col_data[key] = lfile[key_path][:]
                else:
                    col_data[key] = lfile[key_path][:]
        column_order = sorted(col_data.columns,
                            key = lambda x: GC._IMU_CELL_COLUMN_INDEX[x])
        col_data = col_data[column_order]
        return(col_data)

def get_gene_int_md_univ(uuid, keep_missing = True):
    """Gets the gene universal metadata.

    Gets the gene universal metadata from the given dataset
    as a Pandas DataFrame. Sets index preferentially
    as Accession, Gene, or a number index in that order.
    Accession and Gene will not be used if they have non-unique
    values.


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
    with ac__.get_h5_conn(uuid) as lfile:
        row_keys = list(lfile['row_attrs'].keys())
        skip_col = []
        if 'Accession' in row_keys:
            acc = pd.Series(lfile['row_attrs/Accession'])
            if acc.shape == acc.unique().shape:
                row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Accession']))
                skip_col.append('Accession')
            else:
                print('WARNING: internal_metadata.get_gene_int_md_univ()')
                print('  \"Accession\" has non-unique values!')
                row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
            del acc
        elif 'Gene' in row_keys:
            gene = pd.Series(lfile['row_attrs/Gene'])
            if gene.shape == gene.unique().shape:
                row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Gene']))
                skip_col.append('Gene')
            else:
                row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
            del gene
        else:
            print('WARNING: internal_metadata.get_gene_int_md_univ()')
            print('  \"Gene\" and \"Accession\" are both missing!')
            row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
        for key in row_keys:
            if key in skip_col:
                continue
            key_path = 'row_attrs/' + key
            if len(lfile[key_path].shape) == 1:
                if (np.array(list(map(str, list(lfile[key_path])))) == '-1').all():
                    if keep_missing:
                        row_data[key] = lfile[key_path][:]
                else:
                    row_data[key] = lfile[key_path][:]
        # Accession and Gene are handled separately because they
        # aren't handled with the global constants, and so they
        # can't be sorted by their index.
        add_acc = False
        add_gene = False
        if 'Accession' in list(row_data.columns):
            add_acc = True
        if 'Gene' in list(row_data.columns):
            add_gene = True
        columns = list(filter(lambda x: x not in ['Accession', 'Gene'],
                                row_data.columns))
        column_order = sorted(columns,
                                key = lambda x: GC._IMU_GENE_COLUMN_INDEX[x])
        # Prepended in reverse order to ensure that
        # Accession will come before Gene
        if add_gene:
            column_order = ['Gene'] + column_order
        if add_acc:
            column_order = ['Accession'] + column_order
        row_data = row_data[column_order]
        return(row_data)

def get_cell_int_md_author_annot(uuid):
    """Gets the cell author-annotated metadata.

    Gets the cell author-annotated metadata from the given
    dataset as a Pandas DataFrame.

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as cells
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    with ac__.get_h5_conn(uuid) as lfile:
        col_keys = list(lfile['cell_author_annot'].keys())
        col_data = pd.DataFrame(index = np.array(lfile['col_attrs/CellID']))
        for key in col_keys:
            key_path = 'cell_author_annot/' + key
            if len(lfile[key_path].shape) == 1:
                col_data[key] = lfile[key_path][:]
        col_order = lfile['cell_author_annot'].attrs['column_order'].split('|')
        col_data = col_data[col_order]
        return(col_data)

def get_gene_int_md_author_annot(uuid):
    """Gets the gene author-annotated metadata.

    Gets the gene author-annotated metadata from the given
    dataset as a Pandas DataFrame. Sets index preferentially
    as Accession, Gene, or a number index in that order.
    Accession and Gene will not be used if they have non-unique
    values.

    Args:
        uuid: String. UUID of the desired dataset

    Returns:
        Pandas DataFrame with the same number of rows as genes
        in the dataset and where each column is an
        author-annotated internal metadata field.
    
    Raises: None
    """
    with ac__.get_h5_conn(uuid) as lfile:
        row_keys = list(lfile['gene_author_annot'].keys())
        univ_row_keys = list(lfile['row_attrs'].keys())
        if 'Accession' in univ_row_keys:
            acc = pd.Series(lfile['row_attrs/Accession'])
            if acc.shape == acc.unique().shape:
                row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Accession']))
            else:
                print('WARNING: internal_metadata.get_gene_int_md_author_annot()')
                print('  \"Accession\" has non-unique values!')
                row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
            del acc
        elif 'Gene' in univ_row_keys:
            gene = pd.Series(lfile['row_attrs/Gene'])
            if gene.shape == gene.unique().shape:
                row_data = pd.DataFrame(index = np.array(lfile['row_attrs/Gene']))
            else:
                row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
            del gene
        else:
            print('WARNING: internal_metadata.get_gene_int_md_author_annot()')
            print('  \"Gene\" and \"Accession\" are both missing!')
            row_data = pd.DataFrame(index = [i for i in range(lfile['matrix'].shape[0])])
        for key in row_keys:
            key_path = 'gene_author_annot/' + key
            if len(lfile[key_path].shape) == 1:
                row_data[key] = lfile[key_path][:]
        col_order = lfile['gene_author_annot'].attrs['column_order'].split('|')
        row_data = row_data[col_order]
        return(row_data)

def set_cell_int_md_author_annot(uuid, df):
    """Sets the cell author-annotated metadata.

    Sets the cell author-annotated metadata for a given dataset
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
        AssertionError: If df has the wrong number of rows.
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if df.shape[0] != lfile['matrix'].shape[1]:
            raise AssertionError('df has the wrong number of rows!')
        # Keys are the columns
        # If a column was deleted to be overwritten, it is backed
        # up to this dictionary. If column is in keys, then it
        # was written to the loom file. If the value of the
        # column-key is not None, then it is a pandas Series
        # and it has been overwritten.
        columns_written = {} 
        column_order = []
        for col in df.columns:
            try:
                if '/' in col:
                    hdf5_col = re.sub(r'/', r'_', col)
                    print('WARNING: internal_metadata.set_cell_int_md_author_annot()\n'
                          '  \"/\" characters are not allowed in column names!\n'
                         f'  \"{col}\" was changed to \"{hdf5_col}\".\n')
                else:
                    hdf5_col = col
                column_order.append(hdf5_col)
                if hdf5_col in lfile['cell_author_annot'].keys():
                    overwrite_col = gu__.get_yes_or_no(f'Column \"{hdf5_col}\" '
                                                        'already exists!\n'
                                                        'Overwrite? (y/n): ')
                    if not overwrite_col:
                        continue
                    else:
                        columns_written[hdf5_col] = pd.Series(lfile[f'cell_author_annot/{hdf5_col}'])
                        del lfile[f'cell_author_annot/{hdf5_col}']
                if df[col].dtype == object:
                    if h5.__version__ == '2.9.0':
                        dset = lfile.create_dataset(f'cell_author_annot/{hdf5_col}',
                                                    (lfile['matrix'].shape[1], ),
                                                    dtype = h5.special_dtype(vlen = str))
                    elif h5.__version__ == '2.10.0':
                        dset = lfile.create_dataset(f'cell_author_annot/{hdf5_col}',
                                                    (lfile['matrix'].shape[1], ),
                                                    dtype = h5.string_dtype())
                    else:
                        raise AssertionError('Behavior under h5py other '
                                             'than v2.9.0 or v2.10.0 is '
                                             'unknown!')
                    if hdf5_col not in columns_written.keys():
                        columns_written[hdf5_col] = None
                    dset[:] = df[col]
                else:
                    dset = lfile.create_dataset(f'cell_author_annot/{hdf5_col}',
                                                data = df[col])
                    if hdf5_col not in columns_written.keys():
                        columns_written[hdf5_col] = None
            except:
                col_warnings = set() 
                for col_del in columns_written:
                    if columns_written[col_del] is None:
                        del lfile[f'cell_author_annot/{col_del}']
                    elif type(columns_written[col_del]) == pd.core.series.Series:
                        if col_del in lfile['cell_author_annot'].keys():
                            del lfile[f'cell_author_annot/{col_del}']
                        if h5.__version__ == '2.9.0':
                            dset = lfile.create_dataset(f'cell_author_annot/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.special_dtype(vlen = str))
                        elif h5.__version__ == '2.10.0':
                            dset = lfile.create_dataset(f'cell_author_annot/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.string_dtype())
                        else:
                            raise AssertionError('Behavior under h5py other '
                                                 'than v2.9.0 or v2.10.0 is '
                                                 'unknown!')
                        dset[:] = columns_written[col_del]
                        col_warnings.add(col_del)
                    else:
                        print('ERROR: columns_written[\'col_del\'] was not '
                            'None or a pandas Series!')
                            
                raise RuntimeError(f'Error in writing column \"{hdf5_col}\"! '
                                    'Operation aborted. The following '
                                    'columns may have been converted to '
                                    'strings:\n'
                                    '    ' + '\n    '.join(col_warnings))
        lfile['cell_author_annot'].attrs['column_order'] = '|'.join(column_order)

def set_gene_int_md_author_annot(uuid, df):
    """Sets the gene author-annotated metadata.

    Sets the gene author-annotated metadata for a given dataset
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
        AssertionError: If df has the wrong number of rows.
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if df.shape[0] != lfile['matrix'].shape[0]:
            raise AssertionError('df has the wrong number of rows!')
        # Keys are the columns
        # If a column was deleted to be overwritten, it is backed
        # up to this dictionary. If column is in keys, then it
        # was written to the loom file. If the value of the
        # column-key is not None, then it is a pandas Series
        # and it has been overwritten.
        columns_written = {} 
        column_order = []
        for col in df.columns:
            try:
                if '/' in col:
                    hdf5_col = re.sub(r'/', r'_', col)
                    print('WARNING: internal_metadata.set_gene_int_md_author_annot()\n'
                          '  \"/\" characters are not allowed in column names!\n'
                         f'  \"{col}\" was changed to \"{hdf5_col}\".\n')
                else:
                    hdf5_col = col
                column_order.append(hdf5_col)
                if hdf5_col in lfile['gene_author_annot'].keys():
                    overwrite_col = gu__.get_yes_or_no(f'Column \"{hdf5_col}\" '
                                                        'already exists!\n'
                                                        'Overwrite? (y/n): ')
                    if not overwrite_col:
                        continue
                    else:
                        columns_written[hdf5_col] = pd.Series(lfile[f'gene_author_annot/{hdf5_col}'])
                        del lfile[f'cell_author_annot/{hdf5_col}']
                if df[col].dtype == object:
                    if h5.__version__ == '2.9.0':
                        dset = lfile.create_dataset(f'gene_author_annot/{hdf5_col}',
                                                    (lfile['matrix'].shape[0], ),
                                                    dtype = h5.special_dtype(vlen = str))
                    elif h5.__version__ == '2.10.0':
                        dset = lfile.create_dataset(f'gene_author_annot/{hdf5_col}',
                                                    (lfile['matrix'].shape[0], ),
                                                    dtype = h5.string_dtype())
                    else:
                        raise AssertionError('Behavior under h5py other '
                                             'than v2.9.0 or v2.10.0 is '
                                             'unknown!')
                    if hdf5_col not in columns_written.keys():
                        columns_written[hdf5_col] = None
                    dset[:] = df[col]
                else:
                    dset = lfile.create_dataset(f'gene_author_annot/{hdf5_col}',
                                                data = df[col])
                    if hdf5_col not in columns_written.keys():
                        columns_written[hdf5_col] = None
            except:
                col_warnings = set() 
                for col_del in columns_written:
                    if columns_written[col_del] is None:
                        del lfile[f'gene_author_annot/{col_del}']
                    elif type(columns_written[col_del]) == pd.core.series.Series:
                        if col_del in lfile['gene_author_annot'].keys():
                            del lfile[f'gene_author_annot/{col_del}']
                        if h5.__version__ == '2.9.0':
                            dset = lfile.create_dataset(f'gene_author_annot/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.special_dtype(vlen = str))
                        elif h5.__version__ == '2.10.0':
                            dset = lfile.create_dataset(f'gene_author_annot/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.string_dtype())
                        else:
                            raise AssertionError('Behavior under h5py other '
                                                 'than v2.9.0 or v2.10.0 is '
                                                 'unknown!')
                        dset[:] = columns_written[col_del]
                        col_warnings.add(col_del)
                    else:
                        print('ERROR: columns_written[\'col_del\'] was not '
                            'None or a pandas Series!')
                            
                raise RuntimeError(f'Error in writing column \"{hdf5_col}\"! '
                                    'Operation aborted. The following '
                                    'columns may have been converted to '
                                    'strings:\n'
                                    '    ' + '\n    '.join(col_warnings))
        lfile['gene_author_annot'].attrs['column_order'] = '|'.join(column_order)

def set_cell_int_md_univ(uuid, df, batch_key):
    """Sets the cell universal metadata.

    Sets the cell universal metadata for a given dataset
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
        AssertionError: If df has the wrong number of rows, or if
            df has invalid columns. Also if column names have
            an invalid character.
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if df.shape[0] != lfile['matrix'].shape[1]:
            raise AssertionError('df has the wrong number of rows!')
        test_cols = np.array(df.columns)
        if not np.isin(test_cols, np.array(list(GC._IMU_CELL_COLUMN_INDEX.keys()))).all():
            raise AssertionError('df has invalid columns!')
        if np.unique(test_cols).shape[0] != df.columns.shape[0]:
            raise AssertionError('df has non-unique columns!')
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
                if '/' in col:
                    raise AssertionError('Universal internal metadata column '
                                         'names may not have \"/\" characters!')
                if col in lfile['col_attrs'].keys():
                    overwrite_col = gu__.get_yes_or_no(f'Column \"{col}\" '
                                                        'already exists!\n'
                                                        'Overwrite? (y/n): ')
                    if not overwrite_col:
                        continue
                    else:
                        columns_written[col] = pd.Series(lfile[f'col_attrs/{col}'])
                        del lfile[f'col_attrs/{col}']
                if df[col].dtype == object:
                    if h5.__version__ == '2.9.0':
                        dset = lfile.create_dataset(f'col_attrs/{col}',
                                                    (lfile['matrix'].shape[1], ),
                                                    dtype = h5.special_dtype(vlen = str))
                    elif h5.__version__ == '2.10.0':
                        dset = lfile.create_dataset(f'col_attrs/{col}',
                                                    (lfile['matrix'].shape[1], ),
                                                    dtype = h5.string_dtype())
                    else:
                        raise AssertionError('Behavior under h5py other '
                                             'than v2.9.0 or v2.10.0 is '
                                             'unknown!')
                    if col not in columns_written.keys():
                        columns_written[col] = None
                    dset[:] = df[col]
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
                        if h5.__version__ == '2.9.0':
                            dset = lfile.create_dataset(f'col_attrs/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.special_dtype(vlen = str))
                        elif h5.__version__ == '2.10.0':
                            dset = lfile.create_dataset(f'col_attrs/{col_del}',
                                                        (lfile['matrix'].shape[1], ),
                                                        dtype = h5.string_dtype())
                        else:
                            raise AssertionError('Behavior under h5py other '
                                                 'than v2.9.0 or v2.10.0 is '
                                                 'unknown!')
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

def set_gene_int_md_univ(uuid, df):
    """Sets the gene universal metadata.

    Sets the gene universal metadata for a given dataset
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
        AssertionError: If df has the wrong number of rows, or if
            df has invalid columns, if column names have
            an invalid character, or if Accession has non-unique
            values.
        RuntimeError: If anything went wrong during the writing
            process. The method attempts to undo all of the
            edits it made before throwing this error.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if df.shape[0] != lfile['matrix'].shape[0]:
            raise AssertionError('df has the wrong number of rows!')
        test_cols = np.array(df.columns)
        if not np.isin(test_cols, np.array(['Gene', 'Accession']
                                         + list(GC._IMU_GENE_COLUMN_INDEX.keys()))).all():
            raise AssertionError('df has invalid columns!')
        if np.unique(test_cols).shape[0] != df.columns.shape[0]:
            raise AssertionError('df has non-unique columns!')
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
                if '/' in col:
                    raise AssertionError('Universal internal metadata column '
                                         'names may not have \"/\" characters!')
                if col == 'Accession':
                    if df[col].shape != df[col].unique().shape:
                        raise AssertionError('\"Accession\" has non-unique '
                                             'values!')
                if col in lfile['row_attrs'].keys():
                    overwrite_col = gu__.get_yes_or_no(f'Column \"{col}\" '
                                                        'already exists!\n'
                                                        'Overwrite? (y/n): ')
                    if not overwrite_col:
                        continue
                    else:
                        columns_written[col] = pd.Series(lfile[f'row_attrs/{col}'])
                        del lfile[f'row_attrs/{col}']
                if df[col].dtype == object:
                    if h5.__version__ == '2.9.0':
                        dset = lfile.create_dataset(f'row_attrs/{col}',
                                                    (lfile['matrix'].shape[0], ),
                                                    dtype = h5.special_dtype(vlen = str))
                    elif h5.__version__ == '2.10.0':
                        dset = lfile.create_dataset(f'row_attrs/{col}',
                                                    (lfile['matrix'].shape[0], ),
                                                    dtype = h5.string_dtype())
                    else:
                        raise AssertionError('Behavior under h5py other '
                                             'than v2.9.0 or v2.10.0 is '
                                             'unknown!')
                    if col not in columns_written.keys():
                        columns_written[col] = None
                    dset[:] = df[col]
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
                        if h5.__version__ == '2.9.0':
                            dset = lfile.create_dataset(f'row_attrs/{col_del}',
                                                        (lfile['matrix'].shape[0], ),
                                                        dtype = h5.special_dtype(vlen = str))
                        elif h5.__version__ == '2.10.0':
                            dset = lfile.create_dataset(f'row_attrs/{col_del}',
                                                        (lfile['matrix'].shape[0], ),
                                                        dtype = h5.string_dtype())
                        else:
                            raise AssertionError('Behavior under h5py other '
                                                 'than v2.9.0 or v2.10.0 is '
                                                 'unknown!')
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

def get_cell_univ_col_desc(uuid, column):
    """Gets the description attribute for a cell universal column.

    For a given column in the internal cell-specific universal
    metadata, get the "description" attribute, if available.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal cell-specific
            universal metadata column.

    Returns:
        The description if available. Otherwise a warning is printed
        and None is returned.

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid) as lfile:
        if column in lfile['col_attrs'].keys():
            if 'description' in lfile[f'col_attrs/{column}'].attrs.keys():
                desc = lfile[f'col_attrs/{column}'].attrs['description']
                if (desc is None) or (desc == ''):
                    print('No description available!')
                    return
                else:
                    return(desc)
            else:
                print('No description available!')
                return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def get_cell_aa_col_desc(uuid, column):
    """Gets the description attribute for a cell author_annot column.

    For a given column in the internal cell-specific author-annotated
    metadata, get the "description" attribute, if available.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal cell-specific
            author-annotated metadata column.

    Returns:
        The description if available. Otherwise a warning is printed
        and None is returned.

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid) as lfile:
        if column in lfile['cell_author_annot'].keys():
            if 'description' in lfile[f'cell_author_annot/{column}'].attrs.keys():
                desc = lfile[f'cell_author_annot/{column}'].attrs['description']
                if (desc is None) or (desc == ''):
                    print('No description available!')
                    return
                else:
                    return(desc)
            else:
                print('No description available!')
                return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def get_gene_univ_col_desc(uuid, column):
    """Gets the description attribute for a gene universal column.

    For a given column in the internal gene-specific universal
    metadata, get the "description" attribute, if available.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal gene-specific
            universal metadata column.

    Returns:
        The description if available. Otherwise a warning is printed
        and None is returned.

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid) as lfile:
        if column in lfile['row_attrs'].keys():
            if 'description' in lfile[f'row_attrs/{column}'].attrs.keys():
                desc = lfile[f'row_attrs/{column}'].attrs['description']
                if (desc is None) or (desc == ''):
                    print('No description available!')
                    return
                else:
                    return(desc)
            else:
                print('No description available!')
                return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def get_gene_aa_col_desc(uuid, column):
    """Gets the description attribute for a gene author-annotated column.

    For a given column in the internal gene-specific author-annotated
    metadata, get the "description" attribute, if available.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal gene-specific
            author-annotated metadata column.

    Returns:
        The description if available. Otherwise a warning is printed
        and None is returned.

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid) as lfile:
        if column in lfile['gene_author_annot'].keys():
            if 'description' in lfile[f'gene_author_annot/{column}'].attrs.keys():
                desc = lfile[f'gene_author_annot/{column}'].attrs['description']
                if (desc is None) or (desc == ''):
                    print('No description available!')
                    return
                else:
                    return(desc)
            else:
                print('No description available!')
                return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def set_cell_univ_col_desc(uuid, column, desc):
    """Sets the description attribute for a cell universal column.

    For a given column in the internal cell-specific universal
    metadata, set the "description" attribute. Provides a reasonably
    robust interactive interface.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal cell-specific
            universal metadata column.
        desc: String. The description to set.

    Returns: None

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if column in lfile['col_attrs'].keys():
            if 'description' in lfile[f'col_attrs/{column}'].attrs.keys():
                loop = True
                while loop:
                    print(f'A description for {column} already exists.')
                    print('Enter a command:\n'
                          '    OVERWRITE : Overwrite the existing description\n'
                          '    CANCEL    : Cancel the operation\n'
                          '    VIEW      : View the existing description before deciding\n')
                    print()
                    user_input = input('> ')
                    if user_input == 'OVERWRITE':
                        loop = False
                    elif user_input == 'CANCEL':
                        return
                    elif user_input == 'VIEW':
                        loop = True
                        print()
                        print('\"\"\"')
                        print(lfile[f'col_attrs/{column}'].attrs['description'])
                        print('\"\"\"')
                        print()
                    else:
                        loop = True
                        print('ERROR: Please enter a valid command!')
                        print()
            lfile[f'col_attrs/{column}'].attrs['desc'] = desc
            return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def set_cell_aa_col_desc(uuid, column, desc):
    """Sets the description attribute for a cell author-annotated column.

    For a given column in the internal cell-specific author-annotated
    metadata, set the "description" attribute. Provides a reasonably
    robust interactive interface.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal cell-specific
            universal metadata column.
        desc: String. The description to set.

    Returns: None

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if column in lfile['cell_author_annot'].keys():
            if 'description' in lfile[f'cell_author_annot/{column}'].attrs.keys():
                # Ask for overwrite
                loop = True
                while loop:
                    print(f'A description for {column} already exists.')
                    print('Enter a command:\n'
                          '    OVERWRITE : Overwrite the existing description\n'
                          '    CANCEL    : Cancel the operation\n'
                          '    VIEW      : View the existing description before deciding\n')
                    print()
                    user_input = input('> ')
                    if user_input == 'OVERWRITE':
                        loop = False
                    elif user_input == 'CANCEL':
                        return
                    elif user_input == 'VIEW':
                        loop = True
                        print()
                        print('\"\"\"')
                        print(lfile[f'cell_author_annot/{column}'].attrs['description'])
                        print('\"\"\"')
                        print()
                    else:
                        loop = True
                        print('ERROR: Please enter a valid command!')
                        print()
            # Write the description
            lfile[f'cell_author_annot/{column}'].attrs['desc'] = desc
            return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def set_gene_univ_col_desc(uuid, column, desc):
    """Sets the description attribute for a gene universal column.

    For a given column in the internal gene-specific universal
    metadata, set the "description" attribute. Provides a reasonably
    robust interactive interface.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal gene-specific
            universal metadata column.
        desc: String. The description to set.

    Returns: None

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if column in lfile['row_attrs'].keys():
            if 'description' in lfile[f'row_attrs/{column}'].attrs.keys():
                # Ask for overwrite
                loop = True
                while loop:
                    print(f'A description for {column} already exists.')
                    print('Enter a command:\n'
                          '    OVERWRITE : Overwrite the existing description\n'
                          '    CANCEL    : Cancel the operation\n'
                          '    VIEW      : View the existing description before deciding\n')
                    print()
                    user_input = input('> ')
                    if user_input == 'OVERWRITE':
                        loop = False
                    elif user_input == 'CANCEL':
                        return
                    elif user_input == 'VIEW':
                        loop = True
                        print()
                        print('\"\"\"')
                        print(lfile[f'row_attrs/{column}'].attrs['description'])
                        print('\"\"\"')
                        print()
                    else:
                        loop = True
                        print('ERROR: Please enter a valid command!')
                        print()
            # Write the description
            lfile[f'row_attrs/{column}'].attrs['desc'] = desc
            return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def set_gene_aa_col_desc(uuid, column, desc):
    """Sets the description attribute for a gene author-annotated column.

    For a given column in the internal gene-specific author-annotated
    metadata, set the "description" attribute. Provides a reasonably
    robust interactive interface.

    Args:
        uuid: String. The UUID of the desired dataset.
        column: String. The name of the desired internal gene-specific
            author-annotated metadata column.
        desc: String. The description to set.

    Returns: None

    Raises:
        ValueError: If column does not exist.
    """
    with ac__.get_h5_conn(uuid, write = True) as lfile:
        if column in lfile['gene_author_annot'].keys():
            if 'description' in lfile[f'gene_author_annot/{column}'].attrs.keys():
                # Ask for overwrite
                loop = True
                while loop:
                    print(f'A description for {column} already exists.')
                    print('Enter a command:\n'
                          '    OVERWRITE : Overwrite the existing description\n'
                          '    CANCEL    : Cancel the operation\n'
                          '    VIEW      : View the existing description before deciding\n')
                    print()
                    user_input = input('> ')
                    if user_input == 'OVERWRITE':
                        loop = False
                    elif user_input == 'CANCEL':
                        return
                    elif user_input == 'VIEW':
                        loop = True
                        print()
                        print('\"\"\"')
                        print(lfile[f'gene_author_annot/{column}'].attrs['description'])
                        print('\"\"\"')
                        print()
                    else:
                        loop = True
                        print('ERROR: Please enter a valid command!')
                        print()
            # Write the description
            lfile[f'gene_author_annot/{column}'].attrs['desc'] = desc
            return
        else:
            raise ValueError(f'\"{column}\" is not a valid column!')

def main():
    pass

if __name__ == '__main__':
    main()