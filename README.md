# midasmlpy
Python implementation of the midasml approach - providing **estimation** and **prediction** methods for high-dimensional mixed-frequency time-series data

## Authors
* Jonas Striaukas - [jstriaukas](https://github.com/jstriaukas)
* Kris Stern - [krisstern](https://github.com/krisstern)

## About

The *midasmlpy* package implements *estimation* and *prediction* methods for high-dimensional mixed-frequency (MIDAS) time-series and panel data in regression models. 
The regularized MIDAS models are estimated using orthogonal (e.g. Legendre) polynomials and the sparse-group LASSO estimator. 
For more information on the *midasmlpy* approach there are references in the footnotes[^1][^2][^3]. 

The package is equipped with the fast implementation of the sparse-group LASSO estimator by means of proximal block coordinate descent. 
High-dimensional mixed frequency time-series data can also be easily manipulated with functions provided in the package.

## The .f90 source files

These can be found at the [midasmlpy/src](./midasmlpy/src) directory.
Compiled `.so` files can be found at [midasmlpy/compiled](./midasmlpy/compiled) directory.
Please note that we have taken the `.f90` source code for sparse-group LASSO from the repo of the R package `sparsegl` on hosted GitHub at https://github.com/dajmcdon/sparsegl.
We have taken their source code in accordance with their GPL-2.0 license as is without any modification as of November 2nd, 2023 UTC.
What can be found at [midasmlpy/src/sparsegl](./midasmlpy/src/sparsegl) are copies of their `.f90` code.

## Software in other languages

- A Julia implementation of the midasml method is available [here](https://github.com/ababii/Pythia.jl).
- A MATLAB implementation of the midasml method is available [here](https://github.com/jstriaukas/midasml_mat).
- An R implementation of the midasml method is available [here](https://github.com/jstriaukas/midasml).

## Installation

To install the midasmlpy package, download the files and at the directory of files run:
```shell
pip install .
```

## Development

To install the midasmlpy package for development, do the following instead in “editable” mode:
```shell
pip install -e .
```

## Testing

To run the tests, run:
```shell
pytest --pyargs src
```

## Remarks

In case you are running the code on a different platform, you can compile the Fortran code <tt>sglfitF.f90</tt> by using <tt>f2py</tt> which is part of <tt>numpy</tt>. 

[^1]: Babii, A., Ghysels, E., & Striaukas, J. Machine learning time series regressions with an application to nowcasting, (2022) *Journal of Business & Economic Statistics*, Volume 40, Issue 3, 1094-1106. https://doi.org/10.1080/07350015.2021.1899933. 

[^2]: Babii, A., Ghysels, E., & Striaukas, J. High-dimensional Granger causality tests with an application to VIX and news, (2022) *Journal of Financial Econometrics*, nbac023. https://doi.org/10.1093/jjfinec/nbac023.

[^3]: Babii, A., R. Ball, Ghysels, E., & Striaukas, J. Machine learning panel data regressions with heavy-tailed dependent data: Theory and application, (2022) *Journal of Econometrics*. https://doi.org/10.1016/j.jeconom.2022.07.001.
