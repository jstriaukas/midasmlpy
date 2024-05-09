import sys
import os
import setuptools
import numpy
from numpy.distutils.core import Extension
from numpy.distutils.core import setup

# Importing `Extension` and `setup` from `numpy.distutils.core`
# try:
#     from numpy.distutils.core import Extension, setup
# except ImportError:
#     sys.exit("install requires: 'numpy'. use pip or easy_install. \n  $ pip install numpy")

_VERSION = "1.0.0"
f_compile_args = ['-ffree-form']

def read(fname):
    with open(os.path.join(os.path.dirname(__file__), fname)) as _in:
        return _in.read()

def get_lib_dir(dylib):
    import subprocess
    from os.path import realpath, dirname

    p = subprocess.Popen("gfortran -print-file-name={}".format(dylib),
                         stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                         shell=True)
    retcode = p.wait()
    if retcode != 0:
        raise Exception("Failed to find {}".format(dylib))

    libdir = dirname(realpath(p.communicate()[0].strip().decode('ascii')))
    return libdir

if sys.platform == 'darwin':
    GFORTRAN_LIB = get_lib_dir('libgfortran.3.dylib')
    QUADMATH_LIB = get_lib_dir('libquadmath.0.dylib')
    ARGS = ["-Wl,-rpath,{}:{}".format(GFORTRAN_LIB, QUADMATH_LIB)]
    f_compile_args += ARGS
    library_dirs = [GFORTRAN_LIB, QUADMATH_LIB]
else:
    library_dirs = None

midasmlpy_lib = Extension(name='sparsegllog_compiled',
                          sources=['midasmlpy/src/sparseglf90/spmatmul.f90','midasmlpy/src/sparseglf90/log_sgl_subfuns.f90','midasmlpy/src/sparseglf90/sgl_subfuns.f90','midasmlpy/src/sparseglf90/sparsegl.f90','midasmlpy/src/sparseglf90/sparsegllog.f90'],
                          extra_f90_compile_args=f_compile_args,
                          library_dirs=library_dirs)

if __name__ == "__main__":
    setup(name="midasmlpy",
          version=_VERSION,
          description="Python wrapper for midasml",
          long_description=read('README.md'),
          author="Jonas Striaukas and Kris Stern",
          author_email="jonas.striaukas@gmail.com",
          url="https://github.com/jstriaukas/midasmlpy",
          install_requires=[
              "numpy>=1.18.5", 
              "scipy>=1.5.2",
              "scikit-learn>=1.4.2" ],
          python_requires=">=3.8.*",
          setup_requires=["setuptools"],
          ext_modules=[midasmlpy_lib],
          packages=['midasmlpy'],
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Environment :: Console',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3',
              'Programming Language :: Python :: 3.6',
              'Programming Language :: Python :: 3.7',
              'Programming Language :: Python :: 3.8',
              'Programming Language :: Python :: 3 :: Only',
              'Operating System :: OS Independent',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
              'Topic :: Scientific/Engineering'
          ])
