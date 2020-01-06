import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import pandas as pd
import datetime as dt
import geo_access as ga__
import general_utils as gu__

def download_devel_data():
    """Downloads development data.

    Scrapes the single-cell studies database from Valentine Svensson's
    website. It also takes a suitable subset of 10 datasets meeting
    certain conditions:

        * Title contains 'single' or 'Single'
        * Reports a number of cells
        * Date is on or after 2019-01-01
        * Measurement is 'RNA-seq'
        * Data location is a GSE ID

    Datetime-stamps the two files and saves them as tsv files in:
        '../devel_data'
    
    Args: None
    Returns: None
    Raises: None
    """
    path = '../devel_data/'
    url = 'http://www.nxn.se/single-cell-studies/data.tsv'

    data = pd.read_csv(url, sep='\t')

    # Datetime stamp on the files for 2 reasons
    #   1. Never overwrite previous versions of the scraped data
    #   2. We know when the data was scraped
    dt_stamp = gu__.get_timestamp(long = False)
    devel_data = data[(data['Title'].str.contains('[Ss]ingle') == True)
                    & (data['Reported cells total'].map(pd.isna) == False)
                    & (data['Date'] >= 20190101)
                    & (data['Measurement'] == 'RNA-seq')
                    & (data['Data location'].str.contains('^GSE') == True)]
    devel_data = devel_data.sample(n=10,
                                replace=False,
                                axis=0)
    data.to_csv(f'{path}valentine_data_full_{dt_stamp}.tsv',
                sep='\t',
                header=True,
                index=False)
    devel_data.to_csv(f'{path}valentine_data_subset_{dt_stamp}.tsv',
                    sep='\t',
                    header=True,
                    index=False)

def import_devel_data(timestamp, subset = True):
    """Imports a version of development data.

    Imports a specific timestamp of the Valentine Svensson development
    data and returns it as a pandas dataframe.
    
    Args:
        timestamp: String indicating the specific version of the
            development data to import. Should be in the strftime
            format '%y%m%d_%H%M%S'. For example, 2019-03-14 13:03:43
            would be '190314_130343'.
        subset: Boolean indicating whether the 10-entry subset or
            the full database should be imported.
    
    Returns:
        A pandas dataframe containing either the entire database or
        a 10-dataset subset. All columns are objects or floats except
        'Date', which turns into a datetime dtype.

    Raises: None
    """
    path = '../devel_data/'
    file_spec = ''
    if subset:
        file_spec = 'subset'
    else:
        file_spec = 'full'
    timestamp = timestamp + '.tsv'
    file_name = path + '_'.join(['valentine_data', file_spec, timestamp])
    df = pd.read_csv(file_name, sep = '\t', thousands = ',')
    # df['Reported cells total'] is not represented as an integer
    # dtype because it has missing values. One of the unaddressed
    # problems in pandas is that you can't have an integer Series
    # that also has missing values. It has to be a float dtype.
    df['Date'] = pd.to_datetime(df['Date'].astype('str'), format = '%Y%m%d')
    return df

def download_devel_data_soft(path, timestamp, gse_ids = None):
    """Downloads GEO SOFT files for datasets in devel_data.

    Takes a set of GEO Series IDs from the Valentine Svensson
    development data, and downloads their SOFT metadata files
    from GEO.

    Args:
        path: String containing the path to which the files should
            be downloaded.
        timestamp: String containing a date of the strftime() format
            '%Y%m%d_%H%M%S', indicating which version of the
            development data to use.
        gse_ids: If None, will use all GSE IDs in the database. If
            'subset', will use all GSE IDs in that subset. Otherwise,
            must be a list of GSE IDs (e.g. GSE11234) that will be
            used as a match for the GSE IDs in the database.
    
    Returns:
        0 if all executions of ga__.get_series_soft_file()
        returned a 0. 1 if any of them returned a 1.

    Raises: None
    """
    if gse_ids == 'subset':
        gse_ids = import_devel_data(timestamp, subset = True)['Data location']
    elif gse_ids is None:
        df = import_devel_data(timestamp, subset = False)
        gse_ids = df['Data location']
        gse_ids = gse_ids[gse_ids.notnull()]
        gse_ids = gse_ids[gse_ids.str.contains('^GSE') == True]
    else:
        df = import_devel_data(timestamp, subset = False)
        gse_ids = df[df['Data location'].isin(gse_ids)]['Data location']
    exit_flag = 0
    for id in gse_ids:
        exit = ga__.get_series_soft_file(id, path = '../devel_data/')
        if exit == 1:
            exit_flag = 1
    return(exit_flag)
    
def main():
    pass

if __name__ == '__main__':
    main()
