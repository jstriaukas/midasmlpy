# How to compile

## Steps
1. Go to the source file directory
    ````shell
    cd midasmlpy/src/sparseglf90
    ````

2. Compile the f90 code with `f2py`
    ````shell
    python -m numpy.f2py spmatmul.f90 log_sgl_subfuns.f90 sgl_subfuns.f90 sparsegl.f90 sparsegllog.f90 -m sparsegllog_module -h sparsegllog_module_1.pyf

    python -m numpy.f2py --fcompiler=gnu95 -c sparsegllog_module_1.pyf spmatmul.f90 log_sgl_subfuns.f90 sgl_subfuns.f90 sparsegl.f90 sparsegllog.f90
    ````

