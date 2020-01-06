import sys
sys.path.append('/home/nfox/projects/single_cell_database/src')
import pandas as pd
import h5py as h5
import global_constants as gc__
import access as ac__

def get_cell_int_md_univ(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_df = {}
    for k, v in lfile['col_attrs'].items():
        if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[1]):
            new_df[k] = pd.Series(v)
    idx = new_df['CellID']
    del new_df['CellID']
    new_df = pd.DataFrame(new_df)
    new_df.index = idx
    return(new_df)

def get_gene_int_md_univ(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_df = {}
    for k, v in lfile['row_attrs'].items():
        if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[0]):
            new_df[k] = pd.Series(v)
    if 'Gene' in lfile['row_attrs'].keys():
        idx = new_df['Gene']
        del new_df['Gene']
    elif 'Accession' in lfile['row_attrs'].keys():
        idx = new_df['Accession']
        del new_df['Accession']
    else:
        raise AssertionError(f'Loom file {uuid} does not have a '
                              'row_attrs/Gene or a row_attrs/Accession!')
    new_df = pd.DataFrame(new_df)
    new_df.index = idx
    return(new_df)

def get_cell_int_md_author_annot(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_df = {}
    for k, v in lfile['col_attrs/author_annot'].items():
        if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[1]):
            new_df[k] = pd.Series(v)
    idx = pd.Series(lfile['col_attrs/CellID'])
    new_df = pd.DataFrame(new_df)
    new_df.index = idx
    return(new_df)

def get_gene_int_md_author_annot(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_df = {}
    for k, v in lfile['row_attrs/author_annot'].items():
        if (len(v.shape) == 1) and (v.shape[0] == lfile['matrix'].shape[0]):
            new_df[k] = pd.Series(v)
    if 'Gene' in lfile['row_attrs'].keys():
        idx = pd.Series(lfile['row_attrs/Gene'])
    elif 'Accession' in lfile['row_attrs'].keys():
        idx = pd.Series(lfile['row_attrs/Accession'])
    else:
        raise AssertionError(f'Loom file {uuid} does not have a '
                              'row_attrs/Gene or a row_attrs/Accession!')
    new_df = pd.DataFrame(new_df)
    new_df.index = idx
    return(new_df)

def get_cell_int_md(uuid, universal = True):
    if universal:
        return(get_cell_int_md_univ(uuid))
    else:
        return(get_cell_int_md_author_annot(uuid))

def get_gene_int_md(uuid, universal = True):
    if universal:
        return(get_gene_int_md_univ(uuid))
    else:
        return(get_gene_int_md_author_annot(uuid))

def get_cell_int_md_shape_univ(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_dict = {}
    for k, v in lfile['col_attrs'].items():
        new_dict[k] = v.shape
    return(new_dict)

def get_gene_int_md_shape_univ(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_dict = {}
    for k, v in lfile['row_attrs'].items():
        new_dict[k] = v.shape
    return(new_dict)

def get_cell_int_md_shape_author_annot(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_dict = {}
    for k, v in lfile['col_attrs/author_annot'].items():
        new_dict[k] = v.shape
    return(new_dict)

def get_gene_int_md_shape_author_annot(uuid):
    lfile = ac__.get_h5_conn(uuid)
    new_dict = {}
    for k, v in lfile['row_attrs/author_annot'].items():
        new_dict[k] = v.shape
    return(new_dict)

def get_cell_int_md_shape(uuid, universal = True):
    if universal:
        return(get_cell_int_md_shape_univ(uuid))
    else:
        return(get_cell_int_md_shape_author_annot(uuid))

def set_cell_int_md_author_annot(uuid, df):
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[1]:
        raise AssertionError('df has the wrong number of rows!')
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['col_attrs/author_annot'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'col_attrs/author_annot/{col}'])
                    del lfile[f'col_attrs/author_annot/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'col_attrs/author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'col_attrs/author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['col_attrs/author_annot'].keys():
                        del lfile[f'col_attrs/author_annot/{col_del}']
                    dset = lfile.create_dataset(f'col_attrs/author_annot/{col_del}',
                                                (lfile['matrix'].shape[1], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))

def set_gene_int_md_author_annot(uuid, df):
    lfile = ac__.get_h5_conn(uuid, write = True)
    if df.shape[0] != lfile['matrix'].shape[0]:
        raise AssertionError('df has the wrong number of rows!')
    # Keys are the columns
    # If a column was deleted to be overwritten, it is backed
    # up to this dictionary. If column is in keys, then it
    # was written to the loom file. If the value of the
    # column-key is not None, then it is a pandas Series
    # and it has been overwritten.
    columns_written = {} 
    for col in df.columns:
        try:
            if col in lfile['row_attrs/author_annot'].keys():
                skip_col = input(f'Column \"{col}\" already exists!\n'
                                'Overwrite? (y/n): ')[0].lower()
                if skip_col == 'n':
                    continue
                else:
                    columns_written[col] = pd.Series(lfile[f'row_attrs/author_annot/{col}'])
                    del lfile[f'col_attrs/author_annot/{col}']
            if df[col].dtype == object:
                dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                                            (lfile['matrix'].shape[1], ),
                                            dtype = h5.string_dtype())
                dset[:] = df[col]
                if col not in columns_written.keys():
                    columns_written[col] = None
            else:
                dset = lfile.create_dataset(f'row_attrs/author_annot/{col}',
                                            data = df[col])
                if col not in columns_written.keys():
                    columns_written[col] = None
        except:
            col_warnings = set() 
            for col_del in columns_written:
                if columns_written[col_del] is None:
                    del lfile[f'row_attrs/author_annot/{col_del}']
                elif type(columns_written[col_del]) == pd.core.series.Series:
                    if col_del in lfile['row_attrs/author_annot'].keys():
                        del lfile[f'row_attrs/author_annot/{col_del}']
                    dset = lfile.create_dataset(f'row_attrs/author_annot/{col_del}',
                                                (lfile['matrix'].shape[1], ),
                                                dtype = h5.string_dtype())
                    dset[:] = columns_written[col_del]
                    col_warnings.add(col_del)
                else:
                    print('ERROR: columns_written[\'col_del\'] was not '
                          'None or a pandas Series!')
                        
            raise RuntimeError(f'Error in writing column \"{col}\"! '
                                'Operation aborted. The following '
                                'columns may have been converted to '
                                'strings:\n'
                                '    ' + '\n    '.join(col_warnings))



def main():
    pass

if __name__ == '__main__':
    main()