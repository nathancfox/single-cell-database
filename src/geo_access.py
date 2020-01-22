"""Functions to scrape GEO Series data.

Functions that enable easy scraping of GEO Series data and
metadata. Also contains some unused, outdated functions
related to building a new folder and creating log files.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
import sys
sys.path.append('/home/scdb_codebase/single_cell_database/src')
import os
import datetime as dt
import general_utils as gu__

def convert_gse_to_folder(gse_id):
    """Converts GEO Series ID to FTP folder.

    Converts a GEO Series ID into the FTP folder that the files are
    located in, on the NCBI FTP servers. Example ftp URL:
    'ftp://ftp.ncbi.nlm.gov/geo/series/GSEnnn/GSE1/soft/GSE1_family.soft.gz'

    The folders are named by converting the last 3 digits of the
    Series ID to 'nnn'. This Python slice notation returns an empty
    string if the string has length <= 3. Thus, even for GSE1, it
    still returns GSEnnn. This also works for longer GSE IDs such as
    GSE12345 returning GSE12nnn.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series
            ID to convert.
    
    Returns:
        String. The name of the NCBI FTP folder holding the
        information for that GEO Series ID, gse_id. Should be of the
        regex form 'GSE[0-9]*nnn/'
    
    Raises: None
    """
    gse_id = gse_id[3:]
    return 'GSE' + gse_id[:-3] + 'nnn/'

def get_series_soft_url(gse_id):
    """Creates FTP URL for GEO Series soft file.

    Given a GEO Series ID, constructs an FTP URL to download
    the associated SOFT file, holding the metadata for that
    GEO Series and its child GEO Samples.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series
            ID to convert.
    
    Returns:
        String. The full FTP URL to the SOFT file associated with
        the GEO Series indicated by gse_id.

    Raises: None
    """
    ftp_url_base = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/'
    gse_folder = convert_gse_to_folder(gse_id)
    ftp_url_full = (ftp_url_base + gse_folder + gse_id
                  + '/soft/' + gse_id + '_family.soft.gz')
    return(ftp_url_full)

def get_series_suppl_url(gse_id):
    """Creates FTP URL for GEO Series supplementary files.

    Given a GEO Series ID, constructs an FTP URL to download
    the associated supplementary files. The URL will point to
    the folder holding the supplementary files, not to any
    particular file.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series
            ID to convert.
    
    Returns:
        String. The full FTP URL to the supplementary files folder
        associated with the GEO Series indicated by gse_id.

    Raises: None
    """
    ftp_url_base = 'ftp://ftp.ncbi.nlm.nih.gov/geo/series/'
    gse_folder = convert_gse_to_folder(gse_id)
    ftp_url_full = (ftp_url_base + gse_folder + gse_id + '/suppl/')
    return(ftp_url_full)

def get_series_soft_file(gse_id, path):
    """Downloads and unzips a GEO Series SOFT file.

    Given a path and a GEO Series ID, downloads the associated
    SOFT file, and gunzips it in the given path.

    NOTE: Creates and removes a temporary bash script file called:
        'get_series_soft_file_temp.sh'
    in the working directory. This method will fail if there is
    already a file called that, to avoid overwriting anything.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series 
            ID to convert.
        path: String containing the path in which to store the GSE
            SOFT file.

    Returns:
        0 if the expected file exists in the expected place.
        1 otherwise.
    
    Raises:
        Error if this method will overwrite a file with its temporary
        bash script.
    """
    ftp_url = get_series_soft_url(gse_id)
    file_name = gse_id + '_family.soft.gz'
    bash_curl_script = (f'wget -q -P {path} '
                        f'\"{ftp_url}\"')
    bash_gunzip_script = (f'gunzip {path}{file_name}')
    if os.path.exists('get_series_soft_file_temp.sh'):
        sys.exit(('File \"get_series_soft_file_temp.sh\" already exists'
                  ' and will be overwritten by a temp file created by'
                  ' this method. Please rename it to something else'
                  ' to save it.'))
    with open('get_series_soft_file_temp.sh', 'w') as f:
        f.write(bash_curl_script)
        f.write('\n')
        f.write(bash_gunzip_script)
        f.write('\n')
    os.system('bash get_series_soft_file_temp.sh')
    os.remove('get_series_soft_file_temp.sh')
    if os.path.exists(f'{path}{gse_id}_family.soft'):
        return(0)
    else:
        return(1)

def get_series_suppl_files(gse_id, path):
    """Downloads and unzips a GEO Series's supplementary files.

    Given a path and a GEO Series ID, downloads the associated
    supplementary files, and gunzips them in the given path. Unlike
    the similar method, get_series_soft_file(), this method does
    not return an exit code, since it does not know the files it
    should be downloading.

    NOTE: Creates and removes a temporary bash script file called:
        'get_series_suppl_file_temp.sh'
    in the working directory. This method will fail if there is
    already a file called that, to avoid overwriting anything.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series
            ID to convert.
        path: String containing the path in which to store the
            GSE supplementary files.

    Returns: None
    
    Raises:
        Error if this method will overwrite a file with its temporary
        bash script.
    """
    ftp_url = get_series_suppl_url(gse_id)
    bash_curl_script = ('wget -q --no-parent --recursive --level=1 '
                       f'--no-directories -P {path} \"{ftp_url}\"')
    bash_gunzip_script = (f'gunzip {path}*gz')
    if os.path.exists('get_series_suppl_file_temp.sh'):
        raise FileExistsError('File \"get_series_suppl_file_temp.sh\" '
                              'already exists and will be overwritten by '
                              'a temp file created by this method. Please '
                              'rename it to something else to save it.')
    with open('get_series_suppl_file_temp.sh', 'w') as f:
        f.write(bash_curl_script)
        f.write('\n')
        f.write(bash_gunzip_script)
        f.write('\n')
    os.system('bash get_series_suppl_file_temp.sh')
    os.remove('get_series_suppl_file_temp.sh')

def build_new_entry_folder(path):
    """Builds a folder for a new entry in the database.

    Given a path, builds a new folder to store files
    related to a given entry in the database. The folder
    is named with a new UUID.

    Args:
        path: String. The path to the location in which the
            new folder will be created.
    
    Returns:
        A tuple with two members:
            1. The new UUID as a string. Also the name of
               the new directory.
            2. The path to the new directory.
    
    Raises: None
    """
    new_id = gu__.get_uuid()
    os.mkdir(os.path.join(path, new_id))
    return((new_id, os.path.join(path, new_id) + '/'))

def build_new_entry_log(new_id, gse_id, path):
    """Builds a new log file for a new entry in the database.

    Starts a new log file for a new entry in the database. Intended
    to be run after build_new_entry_folder(). Writes basic initial
    information.

    Args:
        new_id: String. The UUID assigned to this entry.
        gse_id: String. The GEO Series ID associated with this entry.
        path: String. The path to the folder associated with
            this entry. The folder should be named something
            that matches new_id.

    Returns: None
    Raises: None
    """
    with open(os.path.join(path, 'log.txt'), 'w') as f:
        f.write(f'UUID         : {new_id}\n')
        f.write(f'GEO ID       : {gse_id}\n')
        datetime_stamp = gu__.get_timestamp()
        f.write(f'Date Created : {datetime_stamp}\n')
        temp = '=' * 80
        f.write(f'{temp}\n')

def get_directory_listing(path):
    """Gets a list of all files in a directory.

    Gets a list of strings denoting all files in a given
    directory. This is non-recursive. All sub-directories
    will end with a '/'.

    Args:
        path: String. The path to the folder to be scanned.

    Returns:
        A sorted list of strings, where each string represents
        a file or sub-directory
    
    Raises: None
    """
    listing = []
    with os.scandir(path) as current_dir:
        for entry in current_dir:
            if entry.is_dir():
                listing.append((entry.name + '/'))
            else:
                listing.append(entry.name)
    listing.sort()
    return(listing)

def log_directory_listing(path, old_listing = None):
    """Logs the state of the directory for a given entry.

    Given a directory, meant to be a directory associated with an
    entry in the database, logs the contents of the directory
    to the log.txt file in that directory. If old_listing is
    given, logs the added or deleted files instead of
    logging the current listing.

    Args:
        path: String. Path to the folder that is to be logged. It
            should also contain the log.txt file.
        old_listing: If None, logs the current contents of the
            directory. If not None, should be a tuple with 2 members:

                1. A set of strings where each string is a file or
                   sub-directory. This should represent a previous
                   state of the directory indicated by path. This
                   is ideally generated by
                   set(get_directory_listing()) because the
                   sub-directories will need to end with '/'.
                2. A string. This is the timestamp that will be
                   logged as the time that this old listing was
                   created. Can be in any format, but will look
                   best in this format: 2019-01-01 13:04:32
    
    Returns: None
    Raises: None
    """
    if old_listing is None:
        listing = get_directory_listing(path)
        with open(os.path.join(path, 'log.txt'), 'a') as f:
            idt = ' ' * 4
            f.write('\n')
            datetime_stamp = gu__.get_timestamp()
            f.write(f'> Directory Listing at {datetime_stamp}\n')
            f.write( '  ----------------------------------------')
            for l in listing:
                f.write(f'{idt}{l}\n')
    else:
        old_listing, old_datetime_stamp = old_listing
        current_listing = set(get_directory_listing(path))
        new_files = list(current_listing - old_listing)
        new_files.sort()
        deleted_files = list(old_listing - current_listing)
        deleted_files.sort()
        with open(os.path.join(path, 'log.txt'), 'a') as f:
            idt = ' ' * 4
            f.write('\n')
            datetime_stamp = gu__.get_timestamp()
            f.write(f'> Change in Directory Listing\n')
            f.write('\n')
            f.write(f'{idt}Old : {old_datetime_stamp}\n')
            f.write(f'{idt}New : {datetime_stamp}\n')
            f.write('\n')
            f.write((f'{idt}NOTE: Old Listing and Timestamp are provided '
                     f'manually.\n{idt}      They may be faked by the code '
                     f'calling this function.\n'))
            f.write('\n')
            if len(new_files) == 0 and len(deleted_files) == 0:
                f.write(f'{idt}No changes!\n')
            elif len(new_files) != 0 and len(deleted_files) == 0:
                f.write(f'{idt}New Files\n{idt}---------\n')
                for entry in new_files:
                    f.write(f'{idt}{idt}{entry}\n')
            elif len(new_files) == 0 and len(deleted_files) != 0:
                f.write(f'{idt}Deleted Files\n{idt}-------------\n')
                for entry in deleted_files:
                    f.write(f'{idt}{idt}{entry}\n')
            else:
                f.write(f'{idt}New Files\n{idt}---------\n')
                for entry in new_files:
                    f.write(f'{idt}{idt}{entry}\n')
                f.write('\n')
                f.write(f'{idt}Deleted Files\n{idt}-------------\n')
                for entry in deleted_files:
                    f.write(f'{idt}{idt}{entry}\n')


def download_series_to_db(gse_id, path, log = True):
    """Downloads a GEO Series to the database.

    Downloads the SOFT file and all associated supplementary
    files for a given GEO Series to a new folder and logs
    the files added.

    Args:
        gse_id: String. The GEO Series ID to scrape.
            e.g. GSE12345
        path: String. The path that the new entry's folder
            will be created in.
        log: Boolean. If True, will log the edits made.
             Otherwise, no log will be created.
    
    Returns:
        new_id: String. The name of the new entry's folder
            and also its UUID.

    Raises: None
    """
    new_id, path = build_new_entry_folder(path)
    if log:
        build_new_entry_log(new_id, gse_id, path)
        old_listing = (set(get_directory_listing(path)),
                    gu__.get_timestamp())
    get_series_soft_file(gse_id, path)
    if log:
        log_directory_listing(path, old_listing = old_listing)
        old_listing = (set(get_directory_listing(path)),
                    gu__.get_timestamp())
    get_series_suppl_files(gse_id, path)
    if log:
        log_directory_listing(path, old_listing = old_listing)
    return(new_id)

def get_sample_char_from_soft(soft_file):
    """Get all sample characteristics from a SOFT file.

    Get a dictionary of dictionaries of all sample characteristics
    for each sample in a SOFT file.

    Args:
        soft_file: String. Filepath to the SOFT file to parse.

    Returns:
        A dictionary of dictionaries. The outer dict has sample
        titles as keys and dicts as values. The inner dicts
        have sample characteristic titles as keys and
        sample characteristic values as values.

        Example:
            {'10X43_1_AAAGACGATCTCGC-1':
                {
                    'strain/background': 'CD-1',
                    'genotype/variation': 'wild type'
                }}

    Raises: None
    """
    with open(soft_file, 'r') as f:
        samples = {}
        new_sample = {}
        new_sample_title = ''
        for line in f:
            if line[:7] == '^SAMPLE':
                if new_sample_title != '':
                    samples[new_sample_title] = new_sample
                    new_sample = {}
                    new_sample_title = ''
                new_sample_title = line.split(' = ')[1].strip(' \n\t')
            elif line[:13] == '!Sample_title':
                new_sample_title = line.split(' = ')[1].strip(' \n\t')
            elif line[:27] == '!Sample_characteristics_ch1':
                line = line.split(' = ')[1]
                line = line.split(': ')
                new_sample[line[0].strip(' \n\t')] = line[1].strip(' \n\t')
        if new_sample_title != '':
            samples[new_sample_title] = new_sample
            new_sample = {}
            new_sample_title = ''
        
    return(samples) 

def main():
    pass

if __name__ == '__main__':
    main()