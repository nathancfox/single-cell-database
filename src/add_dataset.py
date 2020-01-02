import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import utils.create_loom as cl__
import utils.external_metadata as em__
import utils.geo_access as ga__

# Global constants
_PATH = '/home/nfox/projects/single_cell_database/database/'
_GSE_ID = 'GSE100320'

# Global variables
_new_uuid = ''

def download_files():
    global _PATH
    global _new_uuid
    _new_uuid = ga__.download_series_to_db(_GSE_ID,
                                           _PATH,
                                           log = False)

def make_loom_data():
    global _PATH
    global _new_uuid
    expr_matrix, barcodes, features = cl__.get_expr_matrix_from_cellranger(_PATH + _new_uuid + '/')
    cl__.create_loom_file(_PATH + _new_uuid + '/' + 'expr_mat.loom',
                          expr_matrix,
                          barcodes,
                          features)
                        
def add_external_metadata():
    loop_dict = {'y': False, 'n': True}
    loop = True
    while loop:
        new_row = em__.get_new_row_input()
        colnames = em__.get_column_names()
        print()
        max_length = 0
        for col in colnames:
            if len(col) > max_length:
                max_length = len(col)
        for i, col in enumerate(colnames):
            print(f'{i + 1:02}. {col:{max_length}} : {new_row[i]:s}')
        print()
        loop = loop_dict[input(u'\u2191 New Entry. Are you sure (y/n): ').lower()[0]]
    em__.append_row(new_row)

def main():
    pass

if __name__ == '__main__':
    main()
    
