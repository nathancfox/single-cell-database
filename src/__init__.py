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
        del module
    except ImportError:
        import_failures.append(mod)
if len(import_failures) != 0:
    import_failures = ', '.join(import_failures)
    raise ImportError('The following dependencies failed to import: '
                      + import_failures)

# All checks are passed and the import setup is happening
from . import access
from .access import *
