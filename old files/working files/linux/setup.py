from setuptools import Extension,setup
from Cython.Build import cythonize
import numpy as np
import os
os.system("apt-get install libblas-dev liblapack-dev")
setup(
    ext_modules=cythonize(Extension(
	  name="midasml2",
	  sources=[
              'pyx/midasml2.pyx',
              'src/wrapper1.cpp',
              'src/wrapper2.cpp',
              'src/fit_sgl.cpp',
              'src/ld_estim.cpp'
              ],
          language='c++',
	  extra_compile_args=['-w'],
          libraries=['blas','lapack']
	  )
	),
	language='c++',
	include_dirs=[np.get_include(),'include/'],
)
