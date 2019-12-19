"""
README
------

The external metadata file is a text, tab-separated-value file
located at the global variable below: _PATH_TO_DATA.

It currently holds the following columns:

    01. Species
    02. Organ
    03. Number of Cells
    04. Condition
    05. Date Generated
    06. Count Format
    07. UMIs
    08. Spikeins
    09. Technology
    10. DOI
    11. Date Integrated
    12. UUID
    13. File Location
    14. Checksum

Column Guide:

    01. Species
        The species that the cells originated from. Should be a
        two word latin name in standard format. e.g. "Homo sapiens"
    02. Organ
        The organ that the cells originated from. Should be all
        lowercase and as simple and high-level as possible.
        e.g. "brain", "kidney", "blood" 
    03. Number of Cells
        The number of cells in this dataset. Should just be a
        simple integer. e.g. "5402"
    04. Condition
        Whether or not these cells fall under a disease condition
        or an experimental condition. Must be one of these two
        values: {"True", "False"}. Should be True if these are
        assumed to be normal cells from a healthy organisms.
        Should be False otherwise. This could be False if the
        organism is unhealthy, if the cells came from a tumor,
        if the organism or cells were subject to an experimental
        condition other than control.
    05. Date Generated
        The date that these data were created. Often, the date
        of the publication reporting this dataset is used.
        Must be in YYYY-MM-DD format. e.g. "2019-01-01"
    06. Count Format
        The format of the expression values in the actual data.
        Are the read counts raw? Normalized? TPM? Should be
        one of the following:
            {
             'raw', 'cpm', 'tpm', 'rpkm',
             'fpkm', 'other: DESCRIPTION'
            }
        where DESCRIPTION can be anything that adequately describes
        the format of the read counts.
    07. UMIs
        Whether or not the read counts were quantified using UMIs
        or not. Must be one of these two values: {"True", "False"}.
    08. Spikeins
        Whether or not spike-in transcripts were included.
        Must be one of these two values: {"True", "False"}.
    09. Technology
        The protocol used to capture and sequence the reads.
        Should be all lowercase and only alphanumeric.
        e.g. "10xchromiumv2", "smartseq2", "celseq"
    10. DOI
        The DOI for the publication that presented these data.
        Must be in hyperlink format.
        e.g. "https://doi.org/10.1038/s41467-019-13056-x"
    11. Date Integrated
        The date that these data were entered into this database.
        Must be in YYYY-MM-DD format. e.g. "2019-01-01"
    12. UUID
        The UUID associated with this entry in the database.
        Should be generated with the python method: uuid.uuid4()
        Must be in the string form of a UUID:
            12345678-1234-5678-1234-567812345678
        where each of the 32 digits is one hexadecimal digit.
        This can be gotten with str(uuid), if uuid is a 
        python class uuid.UUID
    13. File location
        A full, absolute path to the folder holding the loom
        file and all internal metadata for this dataset.
        Must end with a '/'.
        e.g. "/data/single_cell_database/12345678-1234-5678-1234-567812345678/"
    14. Checksum
        A hash value of the file, to be used for integrity checks.
        Format will depend on the particular hash function used.
        For example, if MD5 is used, the value will be a sequence
        of 32 hexadecimal digits.
        e.g. "9e107d9d372bb6826bd81d3542a419d6"

Column names in the actual file should be all lowercase, all
alphanumeric, and with spaces replaced with underscores.
Columns marked with an asterisk are mandatory.

    01. species*
    02. organ
    03. number_of_cells*
    04. condition
    05. date_generated*
    06. count_format*
    07. umis
    08. spikeins
    09. technology
    10. doi*
    11. date_integrated*
    12. uuid*
    13. file_location*
    14. checksum*

"""
import pandas as pd

# Global Variables
_PATH_TO_DATA = '/home/nfox/projects/single_cell_database/sandbox/external_metadata.tsv'

def get_as_dataframe():
    """Gets external metadata as pandas dataframe."""
    df = pd.read_csv(_PATH_TO_DATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df)

def append_row(new_row):
    """Append a new entry to the external metadata database.

    Columns marked with an asterisk are mandatory.

        01. species*
        02. organ
        03. number_of_cells*
        04. condition
        05. date_generated*
        06. count_format*
        07. umis
        08. spikeins
        09. technology
        10. doi*
        11. date_integrated*
        12. uuid*
        13. file_location*
        14. checksum*

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
    for idx in [0, 2, 4, 5, 9, 10, 11, 12, 13]:
        if new_row[idx] == '':
            raise ValueError('Mandatory column is missing!')
    with open(_PATH_TO_DATA, 'a') as f:
        for i in range(len(new_row)):
            f.write(new_row[i])
            if i != (len(new_row) - 1):
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
    df.to_csv(_PATH_TO_DATA,
              sep = '\t',
              header = True,
              index = False)

def get_shape():
    """Get shape of external metadata."""
    df = pd.read_csv(_PATH_TO_DATA,
                     sep = '\t',
                     header = 0,
                     index_col = None)
    return(df.shape)

def get_column_names():
    """Get list of column names of external metadata."""
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

 