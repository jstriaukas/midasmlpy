# How to compile for Python

cd midasmlpy/src/sparsegl
f2py -c spmatmul.f90 log_sgl_subfuns.f90 sgl_subfuns.f90 sparsegl.f90 sparsegllog.f90 -m sparsegllog_module

