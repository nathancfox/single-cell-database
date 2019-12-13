import pandas as pd
import datetime as dt

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