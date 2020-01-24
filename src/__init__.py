# Assert Python 3
import sys
if sys.version_info[0] != 3:
    raise AssertionError('Python must be Version 3!')

# Dependency checks
import importlib
modules = [
            'datetime',
            'gzip',
            'h5py',
            'loompy',
            'numpy',
            'os',
            'pandas',
            'pprint',
            're',
            'scanpy',
            'scipy',
            'shutil',
            'sys',
            'uuid'
          ]
import_failures = []
for mod in modules:
    try:
        module = importlib.import_module(mod)
        if module.__name__ == 'h5py':
            if module.__version__ not in ['2.9.0', '2.10.0']:
                raise ImportError('h5py must be version 2.9.0 or 2.10.0!')
    except ImportError:
        import_failures.append(mod)
if len(import_failures) != 0:
    import_failures = ', '.join(import_failures)
    raise ImportError('The following dependencies failed to import: '
                      + import_failures)

# Cleanup before import setup
del sys
del importlib
del modules
del import_failures
del mod
del module

# All checks are passed and the import setup is happening
from .access import get_anndata
from .access import get_batch_key
from .access import get_cell_author_annot
from .access import get_cell_ids
from .access import get_cell_univ
from .access import get_column_description
from .access import get_expr_mat_names
from .access import get_expr_mat
from .access import get_extern_md
from .access import get_gene_author_annot
from .access import get_gene_ids
from .access import get_gene_univ
from .access import get_h5_conn
from .access import get_loom_conn
from .access import get_loom_filename
from .access import uuid_to_row