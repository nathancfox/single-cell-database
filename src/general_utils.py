"""Helper functions that are generally helpful.

These are utility helper functions that are not specific to this
database project. These could be useful in other projects and
are generally unrelated to each other.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import datetime as dt
import uuid
import os
import gzip
import shutil
import pandas as pd

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
                    one_per_line = False): 
    """Creates pretty string of a list.

    Creates a pretty string of the items in a list within
    the given width, with the given indent in front of each
    line. Unlike pprint, this method does not quote or
    annotate items in any way. This makes it better for
    printing lists for display, rather than for examination.
    Will produce unexpected behavior if the string version
    of a list element contains a new line character.
    
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

    Args:
        list: The list to convert into a pretty string. Items are
            automatically converted to their string representation,
            so, within reason, the list items do not have to be
            strings.
        width: Integer. The maximum character width of each line,
            including the indent.
        indent: A string to prepend to each line. For example, to
            indent the whole output, use '    ' (4 whitespaces).
        one_per_line: Boolean. If True, all elements will be printed
            on their own lines, without commas, regardless of
            element length.
        
    Returns:
        A single string containing appropriate new lines and indents
        that can be printed. Contains a pretty string version of the
        list.

    Raises: None
    """    
    # Building a string in a loop via concatenation is bad practice.
    # Any operations that add to the end of the pretty string are
    # done by adding substrings to a list, then joining at the end.
    # This operation takes linear time, rather than quadratic.
    # However, string concatenations that take place on a single
    # item, such as adding the ", " are left as string concatenations
    # for clarity, since they do not involve the main pretty string.

    line_counter = 0 
    line_width = 0 
    str_out = []
    str_out.append(indent)
    width -= len(indent) 
    for idx, item in enumerate(list): 
        item = str(item) 
        if idx < len(list) - 1 and (not one_per_line): 
            item += ',' 
        if one_per_line:
            str_out.append(item)
            if idx < len(list) - 1:
                str_out.append('\n' + indent)
            continue
        # Only runs if one_per_line is False
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
    return(str_out)

def convert_dict_to_dataframe(dict_of_dicts):
    """Convert dict of dicts to a Pandas DataFrame.

    Convert a dictionary of dictionaries to a Pandas
    DataFrame, where the keys of the outer dict become
    the Index for the resultant DataFrame. The keys
    of the inner dicts must all be identical and
    become the columns of the resultant DataFrame.

    Example:
        >>> test_dict
        {'Sample_1': {'sex': 'F', 'age': 13},
         'Sample_2': {'sex': 'M', 'age': 14}}
        >>> convert_dict_to_dataframe(test_dict)
                  age sex
        Sample_1   13   F
        Sample_2   14   M
        
    Args:
        dict_of_dicts: Dict. The keys of this dict
            will be the Index of the dataframe. The values
            must be dicts with identical keys. The inner keys
            will be the column names of the dataframe. The
            inner values will be the values of the dataframe.
    
    Returns:
        A Pandas DataFrame holding the converted dict of dicts.

    Raises:
        ValueError: If all inner dicts don't have identical keys.
    """
    first_dict = True
    index = []
    for k, v in dict_of_dicts.items():
        if first_dict:
            ref_keys = pd.Series(sorted(list(v.keys())))
            new_dict = {ref_key:[] for ref_key in ref_keys}
            first_dict = False
        index.append(k)
        test_keys = pd.Series(sorted(list(v.keys())))
        if (test_keys.shape[0] != ref_keys.shape[0]
            or (test_keys != ref_keys).any()):
            raise ValueError('Not all inner dicts have the same keys!')
        for k2, v2 in v.items():
            new_dict[k2].append(v2)
    new_df = pd.DataFrame(new_dict, index = index)
    return(new_df)

def get_yes_or_no(input_string):
    """Get a yes or no response from the user.

    Get a yes or no response from the user in a way
    that handles invalid responses. If this code
    doesn't produce a 'y' or a 'n', then the prompt
    prints an error and loops.

    Code: INPUT_FROM_USER[0].lower()

    Args:
        input_string: String. Prompt fed to input()
    
    Returns:
        Boolean. True if answer was yes. False if no.

    Raises: None
    """
    yn_loop = True
    while yn_loop:
        user_input = input(input_string)
        if user_input is None:
            print('ERROR: Please enter yes or no!')
            continue
        if len(user_input) == 0:
            print('ERROR: Please enter yes or no!')
            continue
        user_input = user_input[0].lower()
        if user_input == 'y':
            yn_loop = False
            return(True)
        elif user_input == 'n':
            yn_loop = False
            return(False)
        else:
            print('ERROR: Please enter yes or no!')
            continue

def main():
    pass

if __name__ == '__main__':
    main()