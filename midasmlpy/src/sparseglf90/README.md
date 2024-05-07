# How to compile for Python
# go to directory
cd midasmlpy/src/sparseglf90
# compile
f2py -c spmatmul.f90 log_sgl_subfuns.f90 sgl_subfuns.f90 sparsegl.f90 sparsegllog.f90 -m sparsegllog_module
# mode to compiled folder
mv sparsegllog_module*.so /Users/m.egelundmuller/Documents/GitHub/midasmlpy/midasmlpy/compiled
/path/to/midasmlpy/compiled

