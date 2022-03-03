# External imports
import numpy as np
import math


def lb(degree, a=0, b=1, jmax=None, X=None):
    """Legendre polynomials shifted to [a,b]

    Parameters
    ----------
    degree : int
        Degree of the polynomial.
    a : int
        Lower shift value (i.e., default - 0).
    b : int
        Upper shift value (i.e., default - 1).
    jmax : int
        Number of high-frequency lags.
    X : ndarray
        Optional evaluation grid vector.

    Returns
    -------
    Psi : ndarray
        Weight matrix with Legendre functions up to ``degree``.
    """
    if jmax is not None:
        X = np.linspace(start=0, stop=1, num=jmax)

    if X is None:
        raise ValueError("X is not provided. Either set X or set jmax.")

    n = len(X)
    P = np.ones([n, degree+2])
    Psi = np.ones([n, degree+1]) / math.sqrt(b-a)

    P[:, 1] = 2 * X / (b-a) - ((b+a) / (b-a))

    if degree > 0:
        for i in range(1, degree+1, 1):
            P[:, i+1] = ((2*i+1)/(i+1)) * P[:, 1] * P[:, i] - i/(i+1) * P[:, i-1]
            Psi[:, i] = np.sqrt((2*i + 1) / (b-a)) @ P[:, i]

    return Psi
