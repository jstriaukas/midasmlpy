# import sys
import os
# from distutils.command.build_ext import build_ext
# from distutils.extension import Extension
from distutils.core import setup

_VERSION = "1.0.0"
f_compile_args = ['-ffree-form']


# def read(fname):
#     with open(os.path.join(os.path.dirname(__file__), fname)) as _in:
#         return _in.read()


# def get_lib_dir(dylib):
#     import subprocess
#     from os.path import realpath, dirname

#     p = subprocess.Popen("gfortran -print-file-name={}".format(dylib),
#                          stdout=subprocess.PIPE, stderr=subprocess.PIPE,
#                          shell=True)
#     retcode = p.wait()
#     if retcode != 0:
#         raise Exception("Failed to find {}".format(dylib))

#     libdir = dirname(realpath(p.communicate()[0].strip().decode('ascii')))
#     return libdir


# if sys.platform == 'darwin':
#     GFORTRAN_LIB = get_lib_dir('libgfortran.3.dylib')
#     QUADMATH_LIB = get_lib_dir('libquadmath.0.dylib')
#     ARGS = ["-Wl,-rpath,{}:{}".format(GFORTRAN_LIB, QUADMATH_LIB)]
#     f_compile_args += ARGS
#     library_dirs = [GFORTRAN_LIB, QUADMATH_LIB]
# else:
#     library_dirs = None


# class f2py_Extension(Extension):

#     def __init__(self, name, sourcedirs):
#         Extension.__init__(self, name, sources=[])
#         self.sourcedirs = [os.path.abspath(sourcedir) for sourcedir in sourcedirs]
#         self.dirs = sourcedirs


# class f2py_Build(build_ext):

#     def run(self):
#         for ext in self.extensions:
#             self.build_extension(ext)

#     def build_extension(self, ext):
#         # compile
#         for ind, to_compile in enumerate(ext.sourcedirs):
#             module_loc = os.path.split(ext.dirs[ind])[0]
#             module_name = os.path.split(to_compile)[1].split('.')[0]
#             os.system('cd %s;f2py -c %s -m %s' % (module_loc, to_compile, module_name))


if __name__ == "__main__":
    setup(name="midasmlpy",
          version=_VERSION,
          description="Python wrapper for midasml",
          long_description=read('README.md'),
          author="Jonas Striaukas, Kris Stern, and Marcus Egelund-Müller",
          author_email="jonas.striaukas@gmail.com",
          url="https://github.com/jstriaukas/midasmlpy",
          install_requires=[
              "numpy>=2.0.0",
              "scipy>=1.14.0",
              "scikit-learn>=1.5.1",
              "pandas>=2.2.2",
              "openpyxl>=3.1.5",
              "Cython>=3.0.10"
          ],
          python_requires=">=3.10.0",
          setup_requires=["setuptools"],
          packages=['midasmlpy'],
          classifiers=[
              'Development Status :: 5 - Production/Stable',
              'Environment :: Console',
              'Programming Language :: Python',
              'Programming Language :: Python :: 3',
              'Operating System :: MacOS',
              'Intended Audience :: Developers',
              'Intended Audience :: Science/Research',
              'License :: OSI Approved :: GNU General Public License v2 (GPLv2)',
              'Topic :: Scientific/Engineering'
          ],
        #   ext_modules=[f2py_Extension('sparsegllog_module',
        #                               ['midasmlpy/src/sparseglf90/log_sgl_subfuns.f90',
        #                                'midasmlpy/src/sparseglf90/spmatmul.f90',
        #                                'midasmlpy/src/sparseglf90/sgl_subfuns.f90',
        #                                'midasmlpy/src/sparseglf90/sparsegl.f90',
        #                                'midasmlpy/src/sparseglf90/sparsegllog.f90',
        #                                'midasmlpy/src/sparseglf90/spmatmul.f90',
        #                                'midasmlpy/src/sparseglf90/sglfitF.f90'])],
        #   cmdclass=dict(build_ext=f2py_Build),
          )
