from distutils.core import setup, Extension
setup(
        name = 'pyTSIL',
        version = '1.0',
        ext_modules = [
            Extension(
                'pyTSIL',
                sources = ['pyTSIL.c'],
                libraries = ['tsil'],
                library_dirs = ['/home/martin/FeynTools/TSIL'],
                include_dirs = ['/home/martin/FeynTools/TSIL']
                )
            ]
        )
