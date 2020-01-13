import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import datetime as dt
import uuid
import os
import gzip
import shutil

def get_timestamp(mode = 'both', long = True):
    """Generates string timestamp.

    Generates a string version of a datestamp, a timestamp,
    or both combined. The three versions are in the strftime()
    formats:
        date: '%y%m%d'
        time: '%H%M%S'
        both: '%y%m%d_%H%M%S'
    
    Args:
        mode: String. Must be one of {'both', 'date', 'time'}
        long: Boolean indicating whether to do a long format
            appropriate for log files, or a short format
            appropriate for filenames.

            Examples
            --------
            long  : 2019-01-01 13:03:24
            short : 190101_130324
    
    Returns:
        String timestamp of the date, time, or both. See above
        for format.

    Raises:
        ValueError: Invalid value passed as mode argument
    """
    now = dt.datetime.now()
    if long:
        date_stamp = now.strftime('%Y-%m-%d')
        time_stamp = now.strftime('%H:%M:%S')
        datetime_stamp = date_stamp + ' ' + time_stamp
    else:
        date_stamp = now.strftime('%y%m%d')
        time_stamp = now.strftime('%H%M%S')
        datetime_stamp = date_stamp + '_' + time_stamp
    if mode == 'both':
        return(datetime_stamp)
    elif mode == 'date':
        return(date_stamp)
    elif mode == 'time':
        return(time_stamp)
    else:
        raise ValueError('mode must be one of {\'both\', \'date\', \'time\'}!')

def get_uuid():
    """Generates string version of a new UUID4."""
    return(str(uuid.uuid4()))

def gunzip(files, remove = True):
    """Gunzips files

    Given a list of file paths, gunzips them, and optionally
    removes the .gz files.

    Args:
        files: Either a string, or a list of strings. Each string
            should be a path to a gzipped file to be gunzipped.
        remove: Boolean. If True, the gzipped files will be
            removed, leaving only the gunzipped files.
    
    Returns: None
    Raises: None
    """
    if isinstance(files, str):
        file = files
        files = []
        files.append(file)
    for f in files:
        path, filename = os.path.split(f)
        if filename[-3:] == '.gz':
            new_filename = filename[:-3]
        else:
            new_filename = filename
        with gzip.open(os.path.join(path, filename), 'rb') as f_in:
            with open(os.path.join(path, new_filename), 'wb') as f_out:
                shutil.copyfileobj(f_in, f_out)
        if remove:
            os.remove(os.path.join(path, filename)) 
    
def pretty_str_list(list, width = 50, indent = '',
                    print_str = False): 
    """Creates pretty string of a list.

    Creates a pretty string of the items in a list within
    the given width, with the given indent in front of each
    line. Unlike pprint, this method does not quote or
    annotate items in any way. This makes it better for
    printing lists for display, rather than for examination.
    
    For example:
        >>> test = ['apple', 'banana', 'cherry', 'date',
                   'elderberry', 'fig', 'grapefruit'] 
        >>> test
        ['apple', 'banana', 'cherry', 'date', 'elderberry', 'fig', 'grapefruit']
        >>> print(pprint.pformat(test, width = 30, indent = 4))
        [    'apple',
             'banana',
             'cherry',
             'date',
             'elderberry',
             'fig',
             'grapefruit']
        >>> print(pretty_str_list(test, width = 30, indent = '    '))
            apple, banana, cherry,
            date, elderberry, fig,
            grapefruit

    Building a string in a loop via concatenation is bad practice.
    Any operations that add to the end of the pretty string are
    done by adding substrings to a list, then joining at the end.
    This operation takes linear time, rather than quadratic.
    However, string concatenations that take place on a single
    item, such as adding the ", " are left as string concatenations
    for clarity, since they do not involve the main pretty string.

    Args:
        list: The list to convert into a pretty string. Items are
            automatically converted to their string representation,
            so, within reason, the list items do not have to be
            strings.
        width: Integer. The maximum character width of each line,
            including the indent.
        indent: A string to preprend to each line. For example, to
            indent the whole output, use '    ' (4 whitespaces).
        print_str: Boolean. If True, the pretty string will be printed,
            and None will be returned. If False, the pretty string
            will be returned and nothing will be printed.
        
    Returns:
        A single string containing appropriate new lines and indents
        that can be printed. Contains a pretty string version of the
        list.

    Raises: None
    """    
    line_counter = 0 
    line_width = 0 
    str_out = []
    str_out.append(indent)
    width -= len(indent) 
    for idx, item in enumerate(list): 
        item = str(item) 
        if idx < len(list) - 1: 
            item += ',' 
        if len(item) >= width: 
            if line_counter != 0: 
                str_out.append('\n' + indent) 
            str_out.append(item)
            str_out.append('\n' + indent)
            line_counter += 2 
            line_width = 0 
        else: 
            # + 1 is for the whitespace after the comma.
            if len(item) + 1 + line_width > width: 
                str_out.append('\n' + indent) 
                str_out.append(item)
                line_counter += 1 
                line_width = len(item) 
            else: 
                if line_counter == 0 and line_width == 0: 
                    pass 
                else: 
                    item = ' ' + item 
                str_out.append(item)
                line_width += len(item) 
    str_out = ''.join(str_out)
    if print_str:
        print(str_out)
    else:
        return(str_out)

def main():
    pass

if __name__ == '__main__':
    main()