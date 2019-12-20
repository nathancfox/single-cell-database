"""
README
------

The external metadata file is a text, tab-separated-value file
located at the global variable below: _PATH_TO_DATA.

Column details can be seen in the global dictionary
below: _COLUMN_DESCRIPTIONS

Mandatory columns are indicated in the global
dictionary below: _COLUMN_MANDATORY

Column names in the actual file should be all lowercase, all
alphanumeric, and with spaces replaced with underscores.
"""
import pandas as pd
import os
import re

# Global Variables
_PATH_TO_DATA = '/home/nfox/projects/single_cell_database/sandbox/test_external_metadata.tsv'
_COLUMN_DESCRIPTIONS = {
    'species': 
        (
            '01. Species\n'
            '    The species that the cells originated from. Should be a\n'
            '    two word latin name in standard format. e.g. \"Homo sapiens\"'
        ),
    'organ':
        (
            '02. Organ\n'
            '    The organ that the cells originated from. Should be all\n'
            '    lowercase and as simple and high-level as possible.\n'
            '    e.g. \"brain\", \"kidney\", \"blood\"'
        ),
    'number_of_cells':
        (
            '03. Number of Cells\n'
            '    The number of cells in this dataset. Should just be a\n'
            '    simple integer. e.g. \"5402\"'
        ),
    'condition':
        (
            '04. Condition\n'
            '    Whether or not these cells fall under a disease condition\n'
            '    or an experimental condition. Must be one of these two\n'
            '    values: {\"True\", \"False\"}. Should be True if these are\n'
            '    assumed to be normal cells from a healthy organisms.\n'
            '    Should be False otherwise. This could be False if the\n'
            '    organism is unhealthy, if the cells came from a tumor,\n'
            '    if the organism or cells were subject to an experimental\n'
            '    condition other than control.'
        ),
    'date_generated':
        (
            '05. Date Generated\n'
            '    The date that these data were created. Often, the date\n'
            '    of the publication reporting this dataset is used.\n'
            '    Must be in YYYY-MM-DD format. e.g. \"2019-01-01\"'
        ),
    'count_format':
        (
            '06. Count Format\n'
            '    The format of the expression values in the actual data.\n'
            '    Are the read counts raw? Normalized? TPM? Should be\n'
            '    one of the following:\n'
            '        {\n'
            '         \'raw\', \'cpm\', \'tpm\', \'rpkm\',\n'
            '         \'fpkm\', \'other: DESCRIPTION\'\n'
            '        }\n'
            '    where DESCRIPTION can be anything that adequately describes\n'
            '    the format of the read counts.'
        ),
    'umis':
        (
            '07. UMIs\n'
            '    Whether or not the read counts were quantified using UMIs\n'
            '    or not. Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'spikeins':
        (
            '08. Spike-ins\n'
            '    Whether or not spike-in transcripts were included.\n'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        ),
    'technology':
        (
            '09. Technology\n'
            '    The protocol used to capture and sequence the reads.\n'
            '    Should be all lowercase and only alphanumeric.\n'
            '    e.g. \"10xchromiumv2\", \"smartseq2\", \"celseq\"'
        ),
    'doi':
        (
            '10. DOI\n'
            '    The DOI for the publication that presented these data.\n'
            '    Must be in hyperlink format.\n'
            '    e.g. \"https://doi.org/10.1038/s41467-019-13056-x\"'
        ),
    'accession':
        (
            '11. Accession\n'
            '    The accession number, if available, for the data. This\n'
            '    can be a GEO Series ID, an ArrayExpress ID, or something\n'
            '    else similar.\n'
            '    e.g. \"GSE108291\", \"E-MTAB-7195\"'
        ),
    'date_integrated':
        (
            '12. Date Integrated\n'
            '    The date that these data were entered into this database.\n'
            '    Must be in YYYY-MM-DD format. e.g. \"2019-01-01\"'
        ),
    'uuid':
        (
            '13. UUID\n'
            '    The UUID associated with this entry in the database.\n'
            '    Should be generated with the python method: uuid.uuid4()\n'
            '    Must be in the string form of a UUID:\n'
            '        12345678-1234-5678-1234-567812345678\n'
            '    where each of the 32 digits is one hexadecimal digit.\n'
            '    This can be gotten with str(uuid), if uuid is a \n'
            '    python class uuid.UUID'
        ),
    'file_location':
        (
            '14. File location\n'
            '    A full, absolute path to the folder holding the loom\n'
            '    file and all internal metadata for this dataset.\n'
            '    Must end with a \'/\'.\n'
            '    e.g. \"/data/single_cell_database/12345678-1234-5678-1234-567812345678/\"'
        ),
    'internal':
        (
            '15. Internal\n'
            '    Whether or not the dataset was scraped from a repository. '
            '    False if it is internally generated data.'
            '    Must be one of these two values: {\"True\", \"False\"}.'
        )
}
_COLUMN_INDEX = {
    'species'        : 0,
    'organ'          : 1,
    'number_of_cells': 2,
    'condition'      : 3,
    'date_generated' : 4,
    'count_format'   : 5,
    'umis'           : 6,
    'spikeins'       : 7,
    'technology'     : 8,
    'doi'            : 9,
    'accession'      : 10,
    'date_integrated': 11,
    'uuid'           : 12,
    'file_location'  : 13,
    'internal'       : 14
}
_COLUMN_MANDATORY = {
    'species'        : True,
    'organ'          : False,
    'number_of_cells': True, 
    'condition'      : False,
    'date_generated' : True,
    'count_format'   : True,
    'umis'           : False,
    'spikeins'       : False,
    'technology'     : False,
    'doi'            : False,
    'accession'      : False,
    'date_integrated': True,
    'uuid'           : True,
    'file_location'  : True,
    'internal'       : True
}

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
    global _COLUMN_MANDATORY
    global _COLUMN_INDEX
    global _PATH_TO_DATA
    n_columns = get_shape()[1]
    if len(new_row) != n_columns:
        raise AssertionError(f'new_row must have length {n_columns}!')
    for k, v in _COLUMN_INDEX.items():
        if (_COLUMN_MANDATORY[k]) and (new_row[v] == ''):
            raise ValueError('Mandatory column is missing!')
    with open(_PATH_TO_DATA, 'a') as f:
        for i in range(len(new_row)):
            f.write(new_row[i])
            if i != (len(new_row) - 1):
                f.write('\t')

def get_new_row_input():
    """Get new entry from user.

    Provides a command line, interactive interface to
    collect the values of a new entry from the user.
    The values are not written directly to the
    external metadata file, rather they are returned
    as a list for external code to parse.

    Args: None
    
    Returns:
        A list of strings, in order of the columns
        in the header of the external metadata file.

    Raises:
        AssertionError: If any element of the returned
                        list did not get edited.
    """
    print()
    print('New Entry to External Metadata')
    print('------------------------------')
    print()
    print('For each column, enter your value, then press Enter.')
    print('If you do not have the information, simply press')
    print('Enter without typing anything. If the column is mandatory,')
    print('indicated with an *, the system will not accept an empty value.')
    print()
    print('To execute one of the following options, enter')
    print('the command in any of the prompts.')
    print('    QUIT    : Abort this new entry')
    print('    RESTART : Start over')
    print('    DESC    : Get a column description')
    print()
    input('Press Enter to continue...')
    print()
    global _COLUMN_DESCRIPTIONS
    global _COLUMN_MANDATORY
    loop = True
    colnames = get_column_names()
    col_idx = 0
    new_vals = ['FILL_ME_IN' for i in range(len(colnames))]
    while loop:
        column_name = _COLUMN_DESCRIPTIONS[colnames[col_idx]]
        column_name = column_name.split('\n')[0][4:]
        if _COLUMN_MANDATORY[colnames[col_idx]]:
            print('Column {:02}: {:s}*'.format(col_idx + 1, column_name))
        else:
            print('Column {:02}: {:s}'.format(col_idx + 1, column_name))
        user_input = input('> ')
        if user_input == 'QUIT':
            print('\nExiting...')
            return
        elif user_input == 'RESTART':
            print('\nStarting over...')
            print()
            print('New Entry to External Metadata')
            print('------------------------------')
            print()
            print('For each column, enter your value, then press Enter.')
            print('If you do not have the information, simply press')
            print('Enter without typing anything. However, if the column')
            print('is mandatory, the system will not accept an empty value.')
            print()
            print('To execute one of the following options, enter')
            print('the corresponding command in any of the prompts.')
            print('    Abort Entry             : QUIT')
            print('    Start Over              : RESTART')
            print('    Get a Column Description: DESC')
            print()
            input('Press Enter to continue...')
            print()
            loop = True
            colnames = get_column_names()
            col_idx = 0
            continue
        elif user_input == 'DESC':
            print()
            print(_COLUMN_DESCRIPTIONS[colnames[col_idx]])
            print()
            print()
            continue
        elif _COLUMN_MANDATORY[colnames[col_idx]] and (user_input) == '':
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
    global _PATH_TO_DATA
    if os.path.exists(_PATH_TO_DATA):
        raise AssertionError(f'{_PATH_TO_DATA} already exists!')
    indices = pd.DataFrame({'keys': list(_COLUMN_INDEX.keys()),
                            'values': list(_COLUMN_INDEX.values())})
    indices = indices.sort_values(by = 'values')
    columns = indices['keys']
    with open(_PATH_TO_DATA, 'w') as f:
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
    global _PATH_TO_DATA
    df.to_csv(_PATH_TO_DATA,
              sep = '\t',
              header = True,
              index = False)

def get_as_dataframe():
    """Gets external metadata as pandas dataframe."""
    global _PATH_TO_DATA
    df = pd.read_csv(_PATH_TO_DATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df)

def get_shape():
    """Get shape of external metadata."""
    global _PATH_TO_DATA
    df = pd.read_csv(_PATH_TO_DATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df.shape)

def get_column_names():
    """Get list of column names of external metadata."""
    global _PATH_TO_DATA
    header = pd.read_csv(_PATH_TO_DATA, sep = '\t', header = None, nrows = 1)
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
        * _PATH_TO_DATA points to an existing regular file
        * _PATH_TO_DATA will open successfully
        * _COLUMN_DESCRIPTIONS, _COLUMN_INDEX, and _COLUMN_MANDATORY
          all have the same length
        * _COLUMN_DESCRIPTIONS, _COLUMN_INDEX, and _COLUMN_MANDATORY
          all have the same keys
        * _PATH_TO_DATA is not an empty file
        * The header of _PATH_TO_DATA has the same number of
          columns as keys in _COLUMN_DESCRIPTIONS
        * The header of _PATH_TO_DATA has the same columns
          as keys in _COLUMN_DESCRIPTIONS
        * The header of _PATH_TO_DATA is in the same order
          as the values in _COLUMN_INDEX
        * The column numbers in the values of _COLUMN_DESCRIPTIONS
          match the values in _COLUMN_INDEX
        * The column titles in the values of _COLUMN_DESCRIPTIONS
          match the keys of _COLUMN_DESCRIPTIONS
    """
    global _PATH_TO_DATA 
    global _COLUMN_DESCRIPTIONS 
    global _COLUMN_INDEX 
    global _COLUMN_MANDATORY 

    if not os.path.exists(_PATH_TO_DATA):
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_DATA does not exist!')
    if not os.path.isfile(_PATH_TO_DATA):
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_DATA is not a '
                             'regular file!')

    if len(_COLUMN_DESCRIPTIONS) != len(_COLUMN_INDEX):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                             '_COLUMN_INDEX are different lengths!')
    if len(_COLUMN_DESCRIPTIONS) != len(_COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_DESCRIPTIONS and '
                             '_COLUMN_MANDATORY are different lengths!')
    if len(_COLUMN_INDEX) != len(_COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _COLUMN_INDEX and '
                             '_COLUMN_MANDATORY are different lengths!')
    desc_keys = list(_COLUMN_DESCRIPTIONS.keys())
    desc_keys.sort()
    indx_keys = list(_COLUMN_INDEX.keys())
    indx_keys.sort()
    mand_keys = list(_COLUMN_MANDATORY.keys())
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
        f = open(_PATH_TO_DATA, 'r')
    except:
        print('GLOBAL VAR ERROR: _PATH_TO_DATA file open failed!')
        return
    header = f.readline()
    f.close()
    if header == '':
        raise AssertionError('GLOBAL VAR ERROR: _PATH_TO_DATA is an empty file!') 
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
    for k, v in _COLUMN_INDEX.items():
        if k != header[v]:
            raise AssertionError('GLOBAL VAR ERROR: The values of '
                                 '_COLUMN_INDEX do not match the '
                                 'order of values in the header '
                                 'at _PATH_TO_VAR!')
        col_desc_index = int(_COLUMN_DESCRIPTIONS[k][:2]) - 1
        if col_desc_index != v:
            raise AssertionError('GLOBAL VAR ERROR: The column numbers '
                                 'in the values of _COLUMN_DESCRIPTIONS '
                                 'do not match the values of _COLUMN_INDEX!')
    for k, v in _COLUMN_DESCRIPTIONS.items():
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
    pass

# Every time this module is run or imported, check
# the integrity of the global constants.
verify_global_constants()

if __name__ == '__main__':
    main()
