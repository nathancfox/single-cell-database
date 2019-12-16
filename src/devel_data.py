import pandas as pd
import datetime as dt

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
    now = dt.datetime.now()
    date_stamp = now.strftime('%y%m%d')
    time_stamp = now.strftime('%H%M%S')

    devel_data = data[(data['Title'].str.contains('[Ss]ingle') == True)
                    & (data['Reported cells total'].map(pd.isna) == False)
                    & (data['Date'] >= 20190101)
                    & (data['Measurement'] == 'RNA-seq')
                    & (data['Data location'].str.contains('^GSE') == True)]
    devel_data = devel_data.sample(n=10,
                                replace=False,
                                axis=0)
    data.to_csv(f'{path}valentine_data_full_{date_stamp}_{time_stamp}.tsv',
                sep='\t',
                header=True,
                index=False)
    devel_data.to_csv(f'{path}valentine_data_subset_{date_stamp}_{time_stamp}.tsv',
                    sep='\t',
                    header=True,
                    index=False)

def import_devel_data(timestamp, subset = True):
    """Imports a version of development data.

    Imports a specific timestamp of the Valentine Svensson development
    data and returns it as a pandas dataframe.
    
    Args:
        timestamp: String indicating the specific version of the development
            data to import. Should be in the strftime format '%y%m%d_%H%M%S'.
            For example, 2019-03-14 13:03:43 would be '190314_130343'.
        subset: Boolean indicating whether the 10-entry subset or the full
            database should be imported.
    
    Returns:
        A pandas dataframe containing either the entire database or
        a 10-dataset subset. All columns are objects except:

            * 'Date' turns into a datetime dtype
            * 'Reported cells total' turns to an int64
    
    Raises:
        Exception: 'File read failed!' on error thrown from pandas.read_csv()
    """
    path = '../devel_data/'
    file_spec = ''
    if subset:
        file_spec = 'subset'
    else:
        file_spec = 'full'
    timestamp = timestamp + '.tsv'
    file_name = path + '_'.join(['valentine_data', file_spec, timestamp])
    try:
        df = pd.read_csv(file_name, sep = '\t')
    except:
        print('File read failed!')
    df['Reported cells total'] = df['Reported cells total'].str.replace(',', '')
    df['Reported cells total'] = df['Reported cells total'].astype('int64')
    df['Date'] = pd.to_datetime(df['Date'].astype('str'), format = '%Y%m%d')
    return df


def main():
    pass

if __name__ == '__main__':
    main()
