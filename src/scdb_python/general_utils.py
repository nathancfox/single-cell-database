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
import tarfile
import pandas as pd
import string

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
        
def extract_tar(tar_file, outpath = "", remove = True):
    with tarfile.open(tar_file) as tar:
        tar_member = tar.next()
        skipped_members = []
        while True:
            if tar_member is None:
                break
            forbidden_strings = ('..', '/')
            if any(forb in tar_member.name for forb in forbidden_strings):
                skipped_members.append(tar_member.name)
                continue
            else:
                tar.extract(tar_member, path = outpath)
            tar_member = tar.next()
    if remove:
        os.remove(tar_file)
    for sm in skipped_members:
        print(f'WARNING: {sm} not extracted due to potentially dangerous '
               'filename. Please extract manually!')
    
def pretty_str_list(list, width = 50, indent = '',
                    sep = ', ', one_per_line = False): 
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
        sep: A string to use to separate entries on a single line.
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
    # sep_end is all trailing whitespace in sep
    # sep_beg is everything before that
    sep_beg = ''
    sep_end = sep
    for i in range(len(sep) - 1, -1, -1):
        if sep[i] in string.whitespace:
            continue
        else:
            sep_beg = sep[:i + 1]
            if i == len(sep) - 1:
                sep_end = ''
            else:
                sep_end = sep[i + 1:]
            break
    for idx, item in enumerate(list): 
        item = str(item) 
        if idx < len(list) - 1 and (not one_per_line): 
            item += sep_beg
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
            if line_width == 0:
                # Because of above if statement, I know that the
                # len(item) is currently < width
                str_out.append(item)
                line_width += len(item)
            elif len(sep_end) + len(item) + line_width > width:
                str_out.append('\n' + indent)
                str_out.append(item)
                line_counter += 1
                line_width = len(item)
            else:
                if line_counter == 0 and line_width == 0:
                    pass
                else:
                    item = sep_end + item
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

def get_yes_or_no(prompt):
    """Get a yes or no response from the user.

    Get a yes or no response from the user in a way
    that handles invalid responses. If this code
    doesn't produce a 'y' or a 'n', then the prompt
    prints an error and loops.

    Code: INPUT_FROM_USER[0].lower()

    Args:
        prompt: String. Prompt fed to input()
    
    Returns:
        Boolean. True if answer was yes. False if no.

    Raises: None
    """
    while True:
        user_input = input(prompt)
        if user_input == '':
            print('ERROR: Please enter yes or no!')
            continue
        if len(user_input) == 0:
            print('ERROR: Please enter yes or no!')
            continue
        user_input = user_input[0].lower()
        if user_input == 'y':
            return(True)
        elif user_input == 'n':
            return(False)
        else:
            print('ERROR: Please enter yes or no!')
            continue

def get_command(options, prompt, first_letter=False,
                case_sensitive=True, error_message=None,
                empty_str_error_msg=None):
    """Get a command from the user in a CLI.

    In a command line interactive application, get a command from
    the user from preset list of options. Provides a stable loop
    that handles bad input.

    Args:
        options: List of Strings. Each string is a valid input option
            to receive from the user. The empty string can be included
            here as a valid option.
        prompt: String. The string used as the prompt for input()
        first_letter: Boolean. If True and all options start with a
            different letter, valid input can be given as the first
            letter of any option. To avoid confusing interfaces,
            'a' and 'A' are not considered different first letters,
            even if the case_sensitive flag is set.
        case_sensitive: Boolean. If True, the input does not have
            to match case to the given options.
        error_message: String. If not None, this will be used as
            a generic error message for all invalid input except
            empty string input. If None, the following will be used:
                'ERROR: Invalid option!'
        empty_str_error_msg: String. If not None, this will be used
            as the error message for invalid empty string input.
            If None, the following will be used:
                'ERROR: Empty string is not a valid option!'
        
    Returns:
        A single string matching one of the options passed in
        options, even if a first_letter abbreviation was
        given by the user or if input is not case-sensitive.

    Raises:
        Assertion Error: If first_letter is set True and there are
            non-unique first letters. This ignores case and the
            case_sensitive flag.
    """
    if error_message is None:
        error_message = 'ERROR: Invalid option!'
    if empty_str_error_msg is None:
        empty_str_error_msg = 'ERROR: Empty string is not a valid option!'
    if first_letter:
        first_letters = set()
        for opt in options:
            if opt == '':
                first_letters.add(opt)
            else:
                first_letters.add(opt[0].lower())
        if len(first_letters) != len(options):
            raise AssertionError('first_letter cannot be True if there are '
                                 'non-unique first letters in the options! '
                                 'Case does not count and there cannot be '
                                 'more than one empty string.')
        del first_letters
    opt_set = set(options)
    if len(opt_set) != len(options):
        raise AssertionError('options has non-unique entries!')
    options = opt_set
    del opt_set
    if '' in options:
        empty_string = True
        options.remove('')
    else:
        empty_string = False
    while True:
        user_input = input(prompt)
        # Take care of empty string case to simplify later code
        if user_input == '':
            if empty_string:
                return('')
            else:
                print(empty_str_error_msg)
                continue
        if user_input in options:
            return(user_input)
        elif not case_sensitive:
            user_input_lower = user_input.lower()
            for opt in options:
                opt_lower = opt.lower()
                if user_input_lower == opt_lower:
                    return(opt)
                elif first_letter and len(user_input) == 1:
                    if user_input_lower == opt_lower[0]:
                        return(opt)
            print(error_message)
            continue
        elif first_letter and len(user_input) == 1:
            for opt in options:
                if user_input == opt[0]:
                    return(opt)
            print(error_message)
            continue
        else:
            print(error_message)
            continue

def pretty_str_text(text, width=80, indent=''):
    """Returns pretty string of string.

    Given a string, returns a prettified, printable
    version of the same string. This is a special case
    of pretty_str_list().

    Args:
        text: String. The string to be converted to a
            pretty printable string.
        width: Integer. The maximum character width of
            the pretty string.
        indent: String. A string to prepend to every line
            of the pretty string.
    
    Returns:
        A printable version of the text string. It conforms
        to the width restriction and adds the indent to the
        beginning of each line.

    Raises: None
    """
    pretty_str = pretty_str_list(text.split(' '), width=width, indent=indent,
                                 sep=' ', one_per_line=False)
    return(pretty_str)

def main():
    pass

if __name__ == '__main__':
    main()