"""Functions to manage the external metadata.

Functions to read, write, and manage the external metadata file
for the database. Does not actually store the location of the
external metadata file. This is stored in global_constants.py.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import pandas as pd
import os
import re
import pprint
import requests
from . import global_constants as GC
from . import general_utils as gu__

def append_row(new_row):
    """Append a new entry to the external metadata database.

    Args:
        new_row: list holding all the entries of the new row
            for the external metadata file. Missing data
            should be included as a '-1'.
    
    Returns: None

    Raises:
        AssertionError: If new_row is the wrong length
        ValueError: If mandatory column is missing
    """
    n_columns = get_shape()[1]
    if len(new_row) != n_columns:
        raise AssertionError(f'new_row must have length {n_columns}!')
    for k, v in GC._EM_COLUMN_INDEX.items():
        if (GC._EM_COLUMN_MANDATORY[k]) and (new_row[v] == ''):
            raise ValueError('Mandatory column is missing!')
    with open(GC.get_PATH_TO_METADATA(), 'a') as f:
        # Temporarily removed because I think it's introducing a bug.
        # f.write('\n')
        for i in range(len(new_row)):
            f.write(new_row[i])
            if i != (len(new_row) - 1):
                f.write('\t')

def get_new_row_input(pre_fill = None):
    """Gets new entry from user.

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
        ValueError: If pre_fill doesn't match expectations.
    """
    intro_text = (
                   '\n'
                   'New Entry to External Metadata\n'
                   '------------------------------\n'
                   '\n'
                   'For each column, enter your value, then press Enter.\n'
                   'If you do not have the information, simply press\n'
                   'Enter without typing anything. If the column is mandatory,\n'
                   'indicated with an *, the system will not accept an empty value.\n'
                   'If a dataset has a characteristic that does not match the\n'
                   'given schema, or if the dataset has multiple values for a\n'
                   'column where multiple values are not allowed, enter \"OTHER\".\n'
                   '\n'
                   'To execute one of the following options, enter\n'
                   'the command in any of the prompts.\n'
                   '    QUIT    : Abort this new entry\n'
                   '    RESTART : Start this new entry over\n'
                   '    DESC    : Get a column description\n'
                   '    LIST    : Get a list of possible values. Note that\n'
                   '              not all columns have this option. Those that\n'
                   '              do, have it mentioned in DESC.\n'
                 )
    print(intro_text)
    input('Press Enter to continue...')
    print()
    loop = True
    colnames = get_column_names()
    col_idx = 0
    new_vals = ['FILL_ME_IN' for i in range(len(colnames))]
    while loop:
        column_name = GC._EM_COLUMN_DESCRIPTIONS[colnames[col_idx]]
        column_name = column_name.split('\n')[0][4:]
        if GC._EM_COLUMN_MANDATORY[colnames[col_idx]]:
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
            print(GC._EM_COLUMN_DESCRIPTIONS[colnames[col_idx]])
            print()
            print()
            continue
        elif user_input == 'LIST':
            if col_idx == 1:
                # tissue
                print()
                for k, v in GC._TISSUE_LIST.items():
                    print(f'{k}')
                    print('-' * (len(k)))
                    print(gu__.pretty_str_list(v, width = 50, indent = '  '))
                    print()
            elif col_idx == 5:
                # count_format
                values = ['raw', 'cpm', 'tpm', 'rpkm', 'fpkm']
                print()
                print('Possible Values')
                print('---------------')
                print(gu__.pretty_str_list(values, width = 50, indent = '  '))
                print()
            elif col_idx == 8:
                # technology
                df = get_as_dataframe()
                values = list(df['technology'].unique())
                values.sort()
                print()
                print('Technologies in Database')
                print('------------------------')
                print(gu__.pretty_str_list(values, width = 50, indent = '  '))
                print()
            else:
                print()
                print('ERROR: This column does not have a \"LIST\" option.')
                print()
            continue
        elif GC._EM_COLUMN_MANDATORY[colnames[col_idx]] and (user_input == '' or user_input == '-1'):
            if pre_fill is not None:
                if colnames[col_idx] in pre_fill.keys():
                    raise ValueError(f'{colnames[col_idx]} is mandatory! '
                                      'It cannot be pre-filled with an '
                                      'empty string or \"-1\" missing value!')
            else:
                print()
                print('ERROR: This column is mandatory, please enter a value!')
                print()
            continue
        else:
            if user_input == '':
                user_input = '-1'
            new_vals[col_idx] = user_input
            col_idx += 1
            if col_idx >= len(colnames):
                loop = False
            print()
    if 'FILL_ME_IN' in new_vals:
        raise AssertionError('At least one column did not get a value!')
    return(new_vals)

def write_header():
    """Writes header of new external metadata file."""
    if os.path.exists(GC.get_PATH_TO_METADATA()):
        raise AssertionError(f'{GC.get_PATH_TO_METADATA()} already exists!')
    indices = pd.DataFrame({'keys': list(GC._EM_COLUMN_INDEX.keys()),
                            'values': list(GC._EM_COLUMN_INDEX.values())})
    indices = indices.sort_values(by = 'values')
    columns = indices['keys']
    with open(GC.get_PATH_TO_METADATA(), 'w') as f:
        for i, k in enumerate(columns):
            f.write(k)
            if i < (len(columns) - 1):
                f.write('\t')

def write_new_file(df):
    """Writes dataframe to file.

    Overwrites the existing external metadata file with
    the given dataframe.

    Args:
        df: pandas dataframe. Data to be overwritten on the
            existing external metadata file.
    
    Returns: None
    Raises: None
    """
    warning = (
        'WARNING! You are about to overwrite the external\n'
        '         metadata file!!! Are you 100% sure you\n'
        '         want to do this??? (yes/no)\n'
    )
    proceed = gu__.get_yes_or_no(f'{warning}\n> ')
    if proceed:
        df.to_csv(GC.get_PATH_TO_METADATA(),
                sep = '\t',
                header = True,
                index = False)
    else:
        pass

def get_as_dataframe():
    """Gets external metadata as pandas dataframe."""
    df = pd.read_csv(GC.get_PATH_TO_METADATA(),
                     sep = '\t',
                     header = 0,
                     index_col = None,
                     parse_dates = ['date_generated', 'date_integrated'])
    return(df)

def get_shape():
    """Gets shape of external metadata."""
    df = pd.read_csv(GC.get_PATH_TO_METADATA(),
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df.shape)

def get_column_names():
    """Gets list of column names of external metadata."""
    header = pd.read_csv(GC.get_PATH_TO_METADATA(), sep = '\t', header = None, nrows = 1)
    header = list(header.iloc[0])
    return(header)

def uuid_to_row(uuid_key, columns = None):
    """Gets values of a row, indexed by UUID.

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
            uuid is a Python uuid.UUID class.
            e.g. "12345678-1234-5678-1234-567812345678"
        columns: list of strings. Each string should be
            a column from the external metadata file
            to include in the returned values. If None,
            the entire row will be returned.
    
    Returns:
        Pandas Series containing the values requested
        that correspond to the uuid_key. Index should
        be the columns variable, if passed. Otherwise,
        will be the columns of the dataframe.

    Raises:
        IndexError: if the uuid_key is not a valid lookup.
    """
    df = get_as_dataframe()
    if uuid_key not in list(df['uuid']):
        raise IndexError('UUID not found in metadata file!')
    row = pd.Series(df[df['uuid'] == uuid_key].iloc[0])
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
                    print('WARNING: At least one column is invalid!')
                warning_flag = True
                values.append('NOT_A_VALID_COLUMN')
        values = pd.Series(values, index = columns)
        return(values)

def get_title_and_author(doi):
    """Gets the title and author for a publication.

    Uses DOI resolution to scrape the title and author for
    a publication.

    Args:
        doi: String. The DOI of the desired publication.

    Returns:
        A 2-member tuple.
            1) A String: the title
            2) A list of strings where each string is an
               author name in the format "GIVEN_NAME FAMILY_NAME".
               If the first author is labeled, it is at the top
               of the list.
    Raises:
        RuntimeError: If the API request returns a non-200
            HTTP status code.
    """
    response = requests.get(f'https://doi.org/{doi}',
                            headers={'Accept': 'application/vnd.citationstyles.csl+json'})
    if response.status_code != 200:
        raise RuntimeError(f'Status Code is {response.status_code}, not 200!')
    response = response.json()
    title = response['title']
    authors = []
    for aut in response['author']:
        if aut['sequence'] == 'first':
            authors.insert(0, aut['given'] + ' ' + aut['family'])
        else:
            authors.append(aut['given'] + ' ' + aut['family'])
    return((title, authors))

def get_abstract(doi):
    """Gets the abstract for a publication in PubMed.

    Uses the NCBI Entrez Utilities to search PubMed for the title
    from a DOI resolution. If there is 1 result, it uses the
    NCBI Entrez Utilities to retrieve the abstract from the PubMed
    entry. Note that this is a bit loose. NCBI's eFetch doesn't
    guarantee an abstract. It returns a long text string that
    Should contain the abstract, but there aren't helpful
    delimiters to guarantee the identity of the scraped text.

    Args:
        doi: String. The DOI of the desired publication. Will be
            resolved to get the title used in the PubMed search.

    Returns:
        A string containing the abstract.

    Raises:
        RuntimeError: If any of the API requests returns a non-200
            HTTP status code.
        AssertionError: If the paragraph before the "abstract"
            violates a light assertion error.
    """
    title = get_title_and_author(doi)[0]
    base_url = 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/'
    esearch_options = (
                          'esearch.fcgi?'
                          'db=pubmed'
                         f'&term={title}[title]'
                          '&retmode=json'
                      )
    esearch_response = requests.get(f'{base_url}{esearch_options}')
    if esearch_response.status_code != 200:
        raise RuntimeError(f'eSearch returned {esearch_response.status_code}'
                            ', not 200!')
    esearch_results = esearch_response.json()
    if str(esearch_results['esearchresult']['count']) != '1':
        raise RuntimeError('eSearch returned '
                          f'{esearch_results["esearchresult"]["count"]} '
                           'results, not 1! You should parse manually.')
    result_id = esearch_results['esearchresult']['idlist'][0]
    efetch_options = (
                          'efetch.fcgi?'
                          'db=pubmed'
                         f'&id={result_id}'
                          '&retmode=text'
                          '&rettype=abstract'
                     )
    efetch_response = requests.get(f'{base_url}{efetch_options}')
    if efetch_response.status_code != 200:
        raise RuntimeError(f'eFetch returned {efetch_response.status_code}'
                            ', not 200!')
    efetch_results = efetch_response.text
    efetch_results = efetch_results.split('\n\n')
    efetch_results = list(map(lambda x: x.strip().replace('\n', ' '),
                              efetch_results))
    if not efetch_results[3].startswith('Author information:'):
        raise AssertionError('Result violates an assumption. You '
                             'should parse manually!')
    return(efetch_results)
    # return(efetch_results[4])

def verify_global_constants():
    """Run assertion tests on global variables.
    
    The following assertions are checked:
        * get_PATH_TO_METADATA() points to an existing regular file
        * get_PATH_TO_METADATA() will open successfully
        * _EM_COLUMN_DESCRIPTIONS, _EM_COLUMN_INDEX, and _EM_COLUMN_MANDATORY
          all have the same length
        * _EM_COLUMN_DESCRIPTIONS, _EM_COLUMN_INDEX, and _EM_COLUMN_MANDATORY
          all have the same keys
        * get_PATH_TO_METADATA() is not an empty file
        * The header of get_PATH_TO_METADATA() has the same number of
          columns as keys in _EM_COLUMN_DESCRIPTIONS
        * The header of get_PATH_TO_METADATA() has the same columns
          as keys in _EM_COLUMN_DESCRIPTIONS
        * The header of get_PATH_TO_METADATA() is in the same order
          as the values in _EM_COLUMN_INDEX
        * The column numbers in the values of _EM_COLUMN_DESCRIPTIONS
          match the values in _EM_COLUMN_INDEX
        * The column titles in the values of _EM_COLUMN_DESCRIPTIONS
          match the keys of _EM_COLUMN_DESCRIPTIONS
    """
    if not os.path.exists(GC.get_PATH_TO_METADATA()):
        raise AssertionError('GLOBAL VAR ERROR: get_PATH_TO_METADATA() does not exist!')
    if not os.path.isfile(GC.get_PATH_TO_METADATA()):
        raise AssertionError('GLOBAL VAR ERROR: get_PATH_TO_METADATA() is not a '
                             'regular file!')

    if len(GC._EM_COLUMN_DESCRIPTIONS) != len(GC._EM_COLUMN_INDEX):
        raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_DESCRIPTIONS and '
                             '_EM_COLUMN_INDEX are different lengths!')
    if len(GC._EM_COLUMN_DESCRIPTIONS) != len(GC._EM_COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_DESCRIPTIONS and '
                             '_EM_COLUMN_MANDATORY are different lengths!')
    if len(GC._EM_COLUMN_INDEX) != len(GC._EM_COLUMN_MANDATORY):
        raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_INDEX and '
                             '_EM_COLUMN_MANDATORY are different lengths!')
    desc_keys = list(GC._EM_COLUMN_DESCRIPTIONS.keys())
    desc_keys.sort()
    indx_keys = list(GC._EM_COLUMN_INDEX.keys())
    indx_keys.sort()
    mand_keys = list(GC._EM_COLUMN_MANDATORY.keys())
    mand_keys.sort()
    for i in range(len(desc_keys)):
        if desc_keys[i] != indx_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_DESCRIPTIONS and '
                                '_EM_COLUMN_INDEX have different keys!')
        if desc_keys[i] != mand_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_DESCRIPTIONS and '
                                '_EM_COLUMN_MANDATORY have different keys!')
        if indx_keys[i] != mand_keys[i]:
            raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_INDEX and '
                                '_EM_COLUMN_MANDATORY have different keys!')

    try:
        f = open(GC.get_PATH_TO_METADATA(), 'r')
    except:
        print('GLOBAL VAR ERROR: get_PATH_TO_METADATA() file open failed!')
        return
    with f as f:
        header = f.readline()
    if header == '':
        raise AssertionError('GLOBAL VAR ERROR: get_PATH_TO_METADATA() is an empty file!') 
    header = header.strip().split(sep = '\t')
    if len(header) != len(desc_keys):
        raise AssertionError('GLOBAL VAR ERROR: _EM_COLUMN_DESCRIPTIONS, '
                             '_EM_COLUMN_INDEX, and _EM_COLUMN_MANDATORY '
                             'have a different length than the header '
                             'at get_PATH_TO_METADATA()!')
    for k in header:
        if k not in desc_keys:
            raise AssertionError('GLOBAL VAR ERROR: The header at '
                                 'get_PATH_TO_METADATA() does not match the '
                                 'keys of _EM_COLUMN_DESCRIPTIONS, '
                                 '_EM_COLUMN_INDEX, and _EM_COLUMN_MANDATORY!')
    for k, v in GC._EM_COLUMN_INDEX.items():
        if k != header[v]:
            raise AssertionError('GLOBAL VAR ERROR: The values of '
                                 '_EM_COLUMN_INDEX do not match the '
                                 'order of values in the header '
                                 'at get_PATH_TO_METADATA()!')
        col_desc_index = int(GC._EM_COLUMN_DESCRIPTIONS[k][:2]) - 1
        if col_desc_index != v:
            raise AssertionError('GLOBAL VAR ERROR: The column numbers '
                                 'in the values of _EM_COLUMN_DESCRIPTIONS '
                                 'do not match the values of _EM_COLUMN_INDEX!')
    for k, v in GC._EM_COLUMN_DESCRIPTIONS.items():
        desc_col = v.split('\n')[0][4:].lower()
        desc_col = re.sub(r'[^a-z0-9_ ]', r'', desc_col)
        desc_col = re.sub(r' ', r'_', desc_col)
        if k != desc_col:
            raise AssertionError('GLOBAL VAR ERROR: The column titles '
                                 'in the values of _EM_COLUMN_DESCRIPTIONS '
                                 'do not match the keys of '
                                 '_EM_COLUMN_DESCRIPTIONS! Note that the '
                                 'column titles are parsed to match '
                                 'the keys, so an imperfect parsing engine '
                                 'may be the problem.')

def main():
    pass

# Every time this module is run or imported, check
# the integrity of the global constants.
# verify_global_constants()

if __name__ == '__main__':
    main()
