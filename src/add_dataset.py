import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import os
import create_loom as cl__
import external_metadata as em__
import geo_access as ga__
import general_utils as gu__
import global_constants as GC

# Global variables
_gse_id = ''
_new_uuid = ''

def print_workflow():
    workflow = (
                 '1. set_new_entry(gse_id)\n'
                 '2. download_files()\n'
                 '3. Rename files appropriately\n'
                 '4. make_loom_data(prefix)\n'
                 '5. add_external_metadata()'
               )
    print()
    print(workflow)
def set_new_entry(new_gse_id):
    global _gse_id
    global _new_uuid
    _gse_id = new_gse_id
    _new_uuid = ''

def download_files():
    global _new_uuid
    _new_uuid = ga__.download_series_to_db(_gse_id,
                                           GC._PATH_TO_DATABASE,
                                           log = False)

def make_loom_data(prefix):
    global _new_uuid
    expr_matrix, barcodes, features = cl__.get_expr_matrix_from_cellranger(GC._PATH_TO_DATABASE + _new_uuid + '/',
                                                                           prefix)
    cl__.create_loom_file(GC._PATH_TO_DATABASE + _new_uuid + '/' + 'expr_mat.loom',
                          expr_matrix,
                          barcodes,
                          features)
                        
def add_external_metadata():
    global _gse_id
    global _new_uuid
    if _new_uuid == '':
        raise AssertionError('Do not run this before download_files()!')
    loop_dict = {'y': False, 'n': True}
    loop = True
    while loop:
        pre_fill = {
                     'accession': _gse_id,
                     'date_integrated': gu__.get_timestamp(mode = 'date'),
                     'uuid': _new_uuid,
                     'file_location': os.path.join(GC._PATH_TO_DATABASE, _new_uuid)
                   }
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
        loop = loop_dict[input(u'\u2191 New Entry. Are you sure (y/n): ').lower()[0]]
    em__.append_row(new_row)

def main():
    pass

if __name__ == '__main__':
    main()
    
