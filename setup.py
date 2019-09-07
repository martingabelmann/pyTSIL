from distutils.core import setup, Extension
from os import environ, path

if 'TSILDIR' not in environ:
    print('Environment variable $TSILDIR not set.\nPoint it to the location where TSIL has been built')
    exit(2)

if not path.isdir(environ['TSILDIR']):
    print('$TSILDIR was not set to a valid dir.')
    exit(2)

if not path.isfile(path.join(environ['TSILDIR'],'libtsil.a')):
    print('TSIL was not built.')
    exit(2)

setup(
        name = 'pyTSIL',
        version = '1.0',
        ext_modules = [
            Extension(
                'pyTSIL',
                sources = ['pyTSIL.c'],
                libraries = ['tsil'],
                library_dirs = [environ['TSILDIR']],
                include_dirs = [environ['TSILDIR']]
                )
            ]
        )
