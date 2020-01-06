"""
README
------

All global constants are stored in the global_constants.py file.

The external metadata file is a text, tab-separated-value file
located at the global constant: _PATH_TO_METADATA.

Column details can be seen in the global constant 
dictionary: _COLUMN_DESCRIPTIONS

Mandatory columns are indicated in the global
constant dictionary: _COLUMN_MANDATORY

Column names in the actual file should be all lowercase, all
alphanumeric, and with spaces replaced with underscores.
"""
import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import pandas as pd
import os
import re
import global_constants as GC

def append_row(new_row):
    """Append a new entry to the external metadata database.

    Args:
        new_row: list holding all the entries of the new row
            for the external metadata file. Missing data
            should be included as a ''.
    
    Returns: None

    Raises:
        AssertionError: If new_row is the wrong length
        ValueError: If mandatory column is missing
    """
    n_columns = get_shape()[1]
    if len(new_row) != n_columns:
        raise AssertionError(f'new_row must have length {n_columns}!')
    for k, v in GC._COLUMN_INDEX.items():
        if (GC._COLUMN_MANDATORY[k]) and (new_row[v] == ''):
            raise ValueError('Mandatory column is missing!')
    with open(GC._PATH_TO_METADATA, 'a') as f:
        # Temporarily removed because I think it's introducing a bug.
        # f.write('\n')
        for i in range(len(new_row)):
            f.write(new_row[i])
            if i != (len(new_row) - 1):
                f.write('\t')

def get_new_row_input(pre_fill = None):
    """Get new entry from user.

    Provides a command line, interactive interface to
    collect the values of a new entry from the user.
    The values are not written directly to the
    external metadata file, rather they are returned
    as a list for external code to parse.

    Args:
        pre_fill: Dict. Each key should be a column name in the
            external metadata file. The values are values to
            be automatically filled in, instead of taking the
            value from the user.
    
    Returns:
        A list of strings, in order of the columns
        in the header of the external metadata file.

    Raises:
        AssertionError: If any element of the returned
                        list did not get edited.
    """
    intro_text = (
                   '\n'
                   'New Entry to External Metadata\n'
                   '------------------------------'
                   '\n'
                   'For each column, enter your value, then press Enter.\n'
                   'If you do not have the information, simply press\n'
                   'Enter without typing anything. If the column is mandatory,\n'
                   'indicated with an *, the system will not accept an empty value.\n'
                   '\n'
                   'To execute one of the following options, enter\n'
                   'the command in any of the prompts.\n'
                   '    QUIT    : Abort this new entry\n'
                   '    RESTART : Start over\n'
                   '    DESC    : Get a column description\n'
                 )
    print(intro_text)
    input('Press Enter to continue...')
    print()
    loop = True
    colnames = get_column_names()
    col_idx = 0
    new_vals = ['FILL_ME_IN' for i in range(len(colnames))]
    while loop:
        column_name = GC._COLUMN_DESCRIPTIONS[colnames[col_idx]]
        column_name = column_name.split('\n')[0][4:]
        if GC._COLUMN_MANDATORY[colnames[col_idx]]:
            print(f'Column {col_idx + 1:02d}: {column_name}*')
        else:
            print(f'Column {col_idx + 1:02d}: {column_name}')
        if pre_fill is not None:
            if type(pre_fill) != dict:
                raise ValueError('pre_fill must be dict!')
            if colnames[col_idx] in pre_fill.keys():
                print('**AUTOMATICALLY FILLED IN**')
                user_input = str(pre_fill[colnames[col_idx]])
                if user_input in ['QUIT', 'RESTART', 'DESC']:
                    raise ValueError('Values in pre_fill cannot be commands '
                                     'in [\'QUIT\', \'RESTART\', \'DESC\']!')
            else:
                user_input = input('> ')
        else:
            user_input = input('> ')
        if user_input == 'QUIT':
            print('\nExiting...')
            return
        elif user_input == 'RESTART':
            print('\nStarting over...')
            print(intro_text)
            loop = True
            colnames = get_column_names()
            col_idx = 0
            continue
        elif user_input == 'DESC':
            print()
            print(GC._COLUMN_DESCRIPTIONS[colnames[col_idx]])
            print()
            print()
            continue
        elif GC._COLUMN_MANDATORY[colnames[col_idx]] and (user_input) == '':
            if pre_fill is not None:
                if colnames[col_idx] in pre_fill.keys():
                    raise ValueError(f'{colnames[col_idx]} is mandatory! '
                                      'It cannot be pre-filled with an '
                                      'empty string!')
            else:
                print()
                print('ERROR: This column is mandatory, please enter a value!')
                print()
            continue
        else:
            new_vals[col_idx] = user_input
            col_idx += 1
            if col_idx >= len(colnames):
                loop = False
            print()
    if 'FILL_ME_IN' in new_vals:
        raise AssertionError('At least one column did not get a value!')
    return(new_vals)

def write_header():
    """Write header of new external metadata file."""
    if os.path.exists(GC._PATH_TO_METADATA):
        raise AssertionError(f'{GC._PATH_TO_METADATA} already exists!')
    indices = pd.DataFrame({'keys': list(GC._COLUMN_INDEX.keys()),
                            'values': list(GC._COLUMN_INDEX.values())})
    indices = indices.sort_values(by = 'values')
    columns = indices['keys']
    with open(GC._PATH_TO_METADATA, 'w') as f:
        for i, k in enumerate(columns):
            f.write(k)
            if i < (len(columns) - 1):
                f.write('\t')

def write_new_file(df):
    """Write dataframe to file.

    Overwrites the existing external metadata file with
    the given dataframe.

    Args:
        df: pandas dataframe. Data to be overwritten on the
            existing external metadata file.
    
    Returns: None
    Raises: None
    """
    df.to_csv(GC._PATH_TO_METADATA,
              sep = '\t',
              header = True,
              index = False)

def get_as_dataframe():
    """Gets external metadata as pandas dataframe."""
    df = pd.read_csv(GC._PATH_TO_METADATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df)

def get_shape():
    """Get shape of external metadata."""
    df = pd.read_csv(GC._PATH_TO_METADATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df.shape)

def get_column_names():
    """Get list of column names of external metadata."""
    header = pd.read_csv(GC._PATH_TO_METADATA, sep = '\t', header = None, nrows = 1)
    header = list(header.iloc[0])
    return(header)

def uuid_to_col(uuid_key, columns = None):
    """Get values of a row, indexed by UUID.

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
    df = get_as_dataframe()
    if uuid_key not in df['uuid']:
        raise IndexError('UUID not found in metadata file!')
    row = pd.Series(df[df['uuid'] == uuid_key, :])
    if columns is None:
        return(row)
    else:
        values = []
        warning_flag = False
        for col in columns:
            if col in df.columns:
                values.append(row[col])
            else:
                if not warning_flag:
                    raise Warning('At least one column is invalid!')
                warning_flag = True
                values.append('NOT_A_VALID_COLUMN')
        values = pd.Series(values, index = columns)
        return(values)

def verify_global_constants():
    """Run assertion tests on global variables.
    
    The following assertions are checked:
        * _PATH_TO_METADATA points to an existing regular file
        * _PATH_TO_METADATA will open successfully
        * _COLUMN_DESCRIPTIONS, _COLUMN_INDEX, and _COLUMN_MANDATORY
          all have the same length
        * _COLUMN_DESCRIPTIONS, _COLUMN_INDEX, and _COLUMN_MANDATORY
          all have the same keys
        * _PATH_TO_METADATA is not an empty file
        * The header of _PATH_TO_METADATA has the same number of
          columns as keys in _COLUMN_DESCRIPTIONS
        * The header of _PATH_TO_METADATA has the same columns
          as keys in _COLUMN_DESCRIPTIONS
        * The header of _PATH_TO_METADATA is in the same order
          as the values in _COLUMN_INDEX
        * The column numbers in the values of _COLUMN_DESCRIPTIONS
          match the values in _COLUMN_INDEX
        * The column titles in the values of _COLUMN_DESCRIPTIONS
          match the keys of _COLUMN_DESCRIPTIONS
    """
    if not os.path.exists(GC._PATH_TO_METADATA):
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_METADATA does not exist!')
    if not os.path.isfile(GC._PATH_TO_METADATA):
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_METADATA is not a '
                             'regular file!')

    if len(GC._COLUMN_DESCRIPTIONS) != len(GC._COLUMN_INDEX):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                             '_COLUMN_INDEX are different lengths!')
    if len(GC._COLUMN_DESCRIPTIONS) != len(GC._COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                             '_COLUMN_MANDATORY are different lengths!')
    if len(GC._COLUMN_INDEX) != len(GC._COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_INDEX and '
                             '_COLUMN_MANDATORY are different lengths!')
    desc_keys = list(GC._COLUMN_DESCRIPTIONS.keys())
    desc_keys.sort()
    indx_keys = list(GC._COLUMN_INDEX.keys())
    indx_keys.sort()
    mand_keys = list(GC._COLUMN_MANDATORY.keys())
    mand_keys.sort()
    for i in range(len(desc_keys)):
        if desc_keys[i] != indx_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                                '_COLUMN_INDEX have different keys!')
        if desc_keys[i] != mand_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                                '_COLUMN_MANDATORY have different keys!')
        if indx_keys[i] != mand_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _COLUMN_INDEX and '
                                '_COLUMN_MANDATORY have different keys!')

    try:
        f = open(GC._PATH_TO_METADATA, 'r')
    except:
        print('GLOBAL VAR ERROR: _PATH_TO_METADATA file open failed!')
        return
    header = f.readline()
    f.close()
    if header == '':
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_METADATA is an empty file!') 
    header = header.strip().split(sep = '\t')
    if len(header) != len(desc_keys):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS, '
                             '_COLUMN_INDEX, and _COLUMN_MANDATORY '
                             'have a different length than the header '
                             'at _PATH_TO_VAR!')
    for k in header:
        if k not in desc_keys:
            raise AssertionError('GLOBAL VAR ERROR: The header at '
                                 '_PATH_TO_VAR does not match the '
                                 'keys of _COLUMN_DESCRIPTIONS, '
                                 '_COLUMN_INDEX, and _COLUMN_MANDATORY!')
    for k, v in GC._COLUMN_INDEX.items():
        if k != header[v]:
            raise AssertionError('GLOBAL VAR ERROR: The values of '
                                 '_COLUMN_INDEX do not match the '
                                 'order of values in the header '
                                 'at _PATH_TO_VAR!')
        col_desc_index = int(GC._COLUMN_DESCRIPTIONS[k][:2]) - 1
        if col_desc_index != v:
            raise AssertionError('GLOBAL VAR ERROR: The column numbers '
                                 'in the values of _COLUMN_DESCRIPTIONS '
                                 'do not match the values of _COLUMN_INDEX!')
    for k, v in GC._COLUMN_DESCRIPTIONS.items():
        desc_col = v.split('\n')[0][4:].lower()
        desc_col = re.sub(r'[^a-z0-9_ ]', r'', desc_col)
        desc_col = re.sub(r' ', r'_', desc_col)
        if k != desc_col:
            raise AssertionError('GLOBAL VAR ERROR: The column titles '
                                 'in the values of _COLUMN_DESCRIPTIONS '
                                 'do not match the keys of '
                                 '_COLUMN_DESCRIPTIONS! Note that the '
                                 'column titles are parsed to match '
                                 'the keys, so an imperfect parsing engine '
                                 'may be the problem.')

def main():
    get_new_row_input()

# Every time this module is run or imported, check
# the integrity of the global constants.
verify_global_constants()

if __name__ == '__main__':
    main()