import os
import sys

def convert_gse_to_folder(gse_id):
    """Convert GEO Series ID to FTP folder.

    Converts a GEO Series ID into the FTP folder that the files are
    located in, on the NCBI FTP servers. Example ftp URL:
    'ftp://ftp.ncbi.nlm.gov/geo/series/GSEnnn/GSE1/soft/GSE1_family.soft.gz'

    The folders are named by converting the last 3 digits of the Series ID
    to 'nnn'. This Python slice notation returns an empty string if the string
    has length <= 3. Thus, even for GSE1, it still returns GSEnnn. This also
    works for longer GSE IDs such as GSE12345 returning GSE12nnn.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series ID to
            convert.
    
    Returns:
        String. The name of the NCBI FTP folder holding the information
        for that GEO Series ID, gse_id. Should be of the regex form
        'GSE[0-9]*nnn/'
    
    Raises: None
    """
    gse_id = gse_id[3:]
    return 'GSE' + gse_id[:-3] + 'nnn/'

def get_series_soft_url(gse_id):
    """Create FTP URL for GEO Series soft file.

    Given a GEO Series ID, constructs an FTP URL to download
    the associated SOFT file, holding the metadata for that
    GEO Series and its child GEO Samples.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series ID to
            convert.
    
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
    """Create FTP URL for GEO Series supplementary files.

    Given a GEO Series ID, constructs an FTP URL to download
    the associated supplementary files. The URL will point to
    the folder holding the supplementary files, not to any
    particular file.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series ID to
            convert.
    
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
    """Download and unzip a GEO Series SOFT file.

    Given a path and a GEO Series ID, downloads the associated
    SOFT file, and gunzips it in the given path.

    NOTE: Creates and removes a temporary bash script file called:
        'get_series_soft_file_temp.sh'
    in the working directory. This method will fail if there is
    already a file called that, to avoid overwriting anything.

    Args:
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series ID to
            convert.
        path: String containing the path in which to store the GSE SOFT file.

    Returns:
        0 if the expected file exists in the expected place. 1 otherwise.
    
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
    """Download and unzip a GEO Series's supplementary files.

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
        gse_id: String of the regex form 'GSE[0-9]+'. The GEO Series ID to
            convert.
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
    bash_gunzip_script = (f'gunzip {path}*')
    if os.path.exists('get_series_suppl_file_temp.sh'):
        sys.exit(('File \"get_series_suppl_file_temp.sh\" already exists'
                  ' and will be overwritten by a temp file created by'
                  ' this method. Please rename it to something else'
                  ' to save it.'))
    with open('get_series_suppl_file_temp.sh', 'w') as f:
        f.write(bash_curl_script)
        f.write('\n')
        f.write(bash_gunzip_script)
        f.write('\n')
    os.system('bash get_series_suppl_file_temp.sh')
    os.remove('get_series_suppl_file_temp.sh')
    

    
def main():
    pass

if __name__ == '__main__':
    main()