"""High-level functions to add a new dataset.

Very high-level functions that consolidate all the helper functions
needed to add a new dataset to the database. Also provides easy
access to some resources like a suggested workflow or internal
metadata schema details.

LICENSE: GNU General Public License v3.0 (see LICENSE file)
"""
# import sys
# sys.path.append('/home/scdb_codebase/single_cell_database/src')
import os
import re
import numpy as np
import pandas as pd
from . import access as ac__
from . import create_loom as cl__
from . import external_metadata as em__
from . import geo_access as ga__
from . import general_utils as gu__
from . import global_constants as GC

# Global variables
_gse_id = ''
_new_uuid = ''

def print_workflow():
    """Prints the workflow to add a new dataset."""
    workflow = (
                 '01. ad__.set_new_entry(gse_id)\n'
                 '02. ad__.download_files()\n'
                 '03. Rename files appropriately\n'
                 '04. Create mat, bar, feat using the appropriate method\n'
                 '05. cl__.create_loom_file()\n'
                 '06. Add layers if applicable\n'
                 '07. Scrape internal author_annot metadata\n'
                 '08. im__.set_cell_int_md_author_annot()\n'
                 '09. im__.set_gene_int_md_author_annot()\n'
                 '10. Construct universal metadata\n'
                 '11. im__.set_cell_int_md_univ()\n'
                 '12. ad__.add_external_metadata()'
               )
    print()
    print(workflow)

def print_tissue_list():
    """Prints the Tissue List."""
    print()
    for k, v in GC._TISSUE_LIST.items():
        print(f'{k}')
        print('-' * (len(k)))
        print(gu__.pretty_str_list(v, width = 50, indent = '  '))
        print()

def print_imu(*args):
    """Prints the internal universal cell-specific metadata schema.

    Provides information about the internal universal cell-specific
    metadata. If called with no arguments, a list of columns is
    printed. If called with arguments, should be a list of columns
    or multiple arguments (one per column name). These columns will
    be printed with full descriptions.

    Args:
        *args: 0 or more positional arguments. Can be column names
            or must be a single list of column names.
    
    Returns:
        If 0 arguments, a list of columns is printed. If 1 argument
        and the argument is a list of column names, those columns
        are printed with their full descriptions. If 1 or more
        arguments and those arguments are column names, those
        columns are printed with their full descriptions.

    Raises: None
    """
    if len(args) == 0:
        column_order = sorted(list(GC._IMU_CELL_COLUMN_INDEX.keys()),
                              key = lambda x: GC._IMU_CELL_COLUMN_INDEX[x])
        print()
        print('Columns')
        print('-------')
        for col in column_order:
            print(f'  {col}')
        print()
    else:
        if len(args) == 1 and type(args[0]) == list:
            args = list(args[0])
        column_order = sorted(list(args),
                              key = lambda x: GC._IMU_CELL_COLUMN_INDEX[x])
        print()
        print('Column Descriptions')
        print('-------------------')
        for col in column_order:
            print(f'{GC._IMU_CELL_COLUMN_DESCRIPTIONS[col]}')
            print()

def setup():
    """Takes care of the first few steps of a new entry.

    Interface for getting started on a new entry from a GEO Series.
    Automatically takes care of the first few steps, and prints
    the remaining steps.

    Args: None
    Returns: None
    Raises: None
    """
    print('Add a New Dataset to the Single-Cell Database')
    print('=============================================')
    print('GEO Series ID')
    new_gse_id = input('> ')
    set_new_entry(new_gse_id)
    global _new_uuid
    print(f'New UUID: {_new_uuid}')
    print('\nDownloading files...')
    download_files()
    print('Done!')
    print()
    print('Remaining Steps')
    print('---------------')
    remaining_steps = ['Rename files appropriately',
                       'Create mat, bar, feat using the appropriate method from create_loom',
                       '>>> cl__.create_loom_file()',
                       'Add layers if applicable',
                       'Scrape internal author-annotated metadata',
                       '>>> im__.set_cell_int_md_author_annot()',
                       '>>> im__.set_gene_int_md_author_annot()',
                       'Construct universal metadata',
                       '>>> im__.set_cell_int_md_univ()',
                       '>>> ad__.add_external_metadata()']
    for i, step in enumerate(remaining_steps):
        print(f'  {i + 1:02d}. {step}')
    print()

def set_new_entry(new_gse_id):
    """Sets global constants and generate a UUID."""
    global _gse_id
    global _new_uuid
    _gse_id = new_gse_id
    _new_uuid = gu__.get_uuid()

def download_files():
    """Downloads GEO files for a new entry."""
    global _new_uuid
    ga__.download_series_to_db(_gse_id,
                               _new_uuid,
                               GC.get_PATH_TO_DATABASE(),
                               log = False)

def add_external_metadata():
    """Gets a new entry's external metadata.

    Robust wrapper around external_metadata.get_new_row_input()
    that provides stable getting of a new entry's external
    metadata and stores it in the external metadata file.

    Args: None
    Returns: None
    Raises: None
    """
    global _gse_id
    global _new_uuid
    if _new_uuid == '':
        raise AssertionError('Do not run this without add_dataset._new_uuid!')
    if _gse_id == '':
        raise AssertionError('Do not run this without add_dataset._gse_id!')
    loop_dict = {'y': False, 'n': True}
    loop = True
    while loop:
        with ac__.get_h5_conn(_new_uuid) as lfile:
            species = ';'.join(np.unique(lfile['col_attrs/species'][:]))
            tissue = ';'.join(np.unique(lfile['col_attrs/tissue'][:]))
            num_cells = lfile['matrix'].shape[1]
        aut_clust = not ac__.get_column_allmissing(_new_uuid, 'cluster',
                                                    var='cell',
                                                    metadata='universal')
        doi = input('Enter the DOI:\n> ')
        title, authors = em__.get_title_and_author(doi)
        title = title.lower()
        title = re.sub(r'[^a-z _-]', r'', title)
        authors = ', '.join(authors)
        pre_fill = {
                     'title': title,
                     'authors': authors,
                     'species': species,
                     'tissue': tissue,
                     'number_of_cells': num_cells,
                     'author_clusters': aut_clust,
                     'accession': _gse_id,
                     'date_integrated': gu__.get_timestamp(mode = 'date'),
                     'uuid': _new_uuid,
                   }
        try:
            abstract = em__.get_abstract(doi)
        except:
            abstract = None
        if abstract is not None:
            print('I attempted to scrape this abstract.')
            print('\nAbstract:\n')
            print(gu__.pretty_str_text(abstract, width = 60,
                                       indent='    '))
            print()
            correct = gu__.get_yes_or_no('Is this correct? (y/n)\n> ')
            if correct:
                pre_fill['abstract'] = abstract
            else:
                print('\nMy bad. You can fill it in manually later.\n')
        new_row = em__.get_new_row_input(pre_fill = pre_fill)
        colnames = em__.get_column_names()
        print()
        max_length = 0
        for col in colnames:
            if len(col) > max_length:
                max_length = len(col)
        for i, col in enumerate(colnames):
            print(f'{i + 1:02}. {col:{max_length}} : {new_row[i]:s}')
        print()
        # Unicode-2191 is an up arrow
        if gu__.get_yes_or_no(u'\u2191 New Entry. Are you sure? (y/n): '):
            loop = False
    em__.append_row(new_row)

def add_note(uuid=None):
    """Interactively adds a note about the new entry.

    Args:
        uuid: String. Indicates the UUID the note is for. If None,
            global _new_uuid will be used.

    Returns: None
    Raises: None
    """
    global _new_uuid
    if uuid is None:
        try:
            ac__.get_loom_filename(_new_uuid)
        except ValueError as e:
            raise ValueError('add_dataset._new_uuid is an invalid UUID!') from e
        uuid = _new_uuid
    else:
        try:
            ac__.get_loom_filename(uuid)
        except ValueError as e:
            raise ValueError('uuid is an invalid UUID!') from e
    notes = pd.read_csv(os.path.join(GC.get_PATH_TO_DATABASE(), 'notes.tsv'),
                        sep='\t', header=None, names=['uuid', 'desc'])
    notes.index = notes['uuid']
    if uuid not in list(notes.index):
        print(f'UUID: {uuid}')
        print(f'Enter your note on a single line')
        new_note = input('> ')
        notes = notes.append(pd.DataFrame({'uuid': [uuid], 'desc': [new_note]}), ignore_index=True)
        print(notes)
        notes.to_csv(os.path.join(GC.get_PATH_TO_DATABASE(), 'notes.tsv'),
                     sep='\t', header=False, index=False)
        return
    else:
        print(f'UUID: {uuid}')
        print('This UUID already has a note!')
        while True:
            print('Enter a command.')
            print('    [V]IEW      : View the current note')
            print('    [O]VERWRITE : Overwrite the current note')
            print('    [A]PPEND    : Append to the current note')
            print('    [C]ANCEL    : Cancel this operation')
            print()
            choice = gu__.get_command(['VIEW', 'OVERWRITE', 'APPEND', 'CANCEL'],
                                    prompt='> ', first_letter=True,
                                    case_sensitive=False)
            if choice == 'VIEW':
                print()
                print(gu__.pretty_str_text(notes.loc[uuid, 'desc'],
                                        width=50))
                print()
                continue
            elif choice == 'OVERWRITE':
                print()
                print('Enter your new description on a single line.')
                new_desc = input('> ')
                notes.loc[uuid, 'desc'] = new_desc
                notes.to_csv(os.path.join(GC.get_PATH_TO_DATABASE(), 'notes.tsv'),
                             sep='\t', header=False, index=False)
                return
            elif choice == 'APPEND':
                print()
                print('Enter your addition on a single line.')
                new_desc = input('> ')
                notes.loc[uuid, 'desc'] = notes.loc[uuid, 'desc'] + new_desc
                notes.to_csv(os.path.join(GC.get_PATH_TO_DATABASE(), 'notes.tsv'),
                             sep='\t', header=False, index=False)
                return
            elif choice == 'CANCEL':
                return
            else:
                raise AssertionError('general_utils.get_command() returned an invalid value!')

def main():
    pass

if __name__ == '__main__':
    main()
    
