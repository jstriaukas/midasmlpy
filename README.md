# midasmlpy
Python implementation of the midasml approach - providing **estimation** and **prediction** methods for high-dimensional mixed-frequency time-series data

## Authors
* Jonas Striaukas - [jstriaukas](https://github.com/jstriaukas)
* Kris Stern - [slim-patchy](https://github.com/slim-patchy)

## About

The *midasmlpy* package implements *estimation* and *prediction* methods for high-dimensional mixed-frequency (MIDAS) time-series and panel data in regression models. 
The regularized MIDAS models are estimated using orthogonal (e.g. Legendre) polynomials and the sparse-group LASSO estimator. 
For more information on the *midasmlpy* approach there are references in the footnotes[^1]. 

The package is equipped with the fast implementation of the sparse-group LASSO estimator by means of proximal block coordinate descent. 
High-dimensional mixed frequency time-series data can also be easily manipulated with functions provided in the package.

## Software in other languages

- Julia implmentation of the midasml method is available [here](https://github.com/ababii/Pythia.jl).
- MATLAB implmentation of the midasml method is available [here](https://github.com/jstriaukas/midasml_mat).
- R implementation of the midasml method is available [here](https://github.com/jstriaukas/midasml).

## Details on the Fortran code

The main subroutines are written in Fortran 90. The code is compiled for several OS: 

- sglfitF.cpython-38-x86_64-linux-gnu.so - compiled for Linux. Version: 20.04.4. Compiled with f2py3 and gfortran compiler.
- sglfitF.cp37-win_amd64.pyd - compiled for Windows. Version: Windows 10. Compiled with f2py and gfortran compiler.

In case you are running the code on a different platform, you can compile the Fortran code  <tt>sglfitF.f90</tt> by using  <tt>f2py</tt> which is part of  <tt>numpy</tt>. 

[^1]: Babii, A., Ghysels, E., & Striaukas, J. (2021). Machine learning time series regressions with an application to nowcasting. 
Journal of Business & Economic Statistics, 1-23.
Now available at https://doi.org/10.1080/07350015.2021.1899933.
