import os
from distutils.core import setup
from distutils.extension import Extension
from Cython.Build import cythonize

mdgxhome = os.environ.get("MDGXHOME")

if mdgxhome is not None:
    library_dirs = [os.path.join(mdgxhome, 'lib'),]
    include_dirs = [os.path.join(mdgxhome, 'include'),]
else:
    amberhome = os.environ.get('AMBERHOME')
    library_dirs = [os.path.join(amberhome, 'lib'),]
    include_dirs = [os.path.join(amberhome, 'include'),]

extra_compile_args = ['-O0', '-ggdb']
extra_link_args = ['-O0', '-ggdb']

cython_directives = {
    'embedsignature': True,
    'boundscheck': False,
    'wraparound': False,
}

ext = Extension(
        "mdgx.mdgx",
        ["mdgx/mdgx.pyx",],
        libraries=['mdgx', 'netcdf',],
        library_dirs=["/home/steboss/Software/amber16/lib"],
        include_dirs=["/home/steboss/Software/amber16/include", "/home/steboss/Simulation/mdgx_python_api/pymdgx/2_readCommFile/mdgx/Vectors",\
	              "/home/steboss/Simulation/mdgx_python_api/pymdgx/2_readCommFile/mdgx/MD",\
                      "/home/steboss/Simulation/mdgx_python_api/pymdgx/2_readCommFile/mdgx/IO"],
        extra_compile_args=extra_compile_args,
        extra_link_args=extra_link_args,
        )

setup(
    name='mdgx',
    packages=['mdgx'],
    version='0.0.1',
    ext_modules=cythonize([ext, ],
        compiler_directives=cython_directives,
    ),
)
