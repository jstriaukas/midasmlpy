# External imports
import numpy as np
import math


def lb(degree, a=0, b=1, jmax=None, X=None):
    """Legendre polynomials shifted to [a,b]

    :param degree:
        Degree of the polynomial.
    :type degree: int
    :param a:
        Lower shift value (i.e., default - 0).
    :type a: int
    :param b:
        Upper shift value (i.e., default - 1).
    :type b: int
    :param jmax:
        Number of high-frequency lags.
    :type jmax: int
    :param X:
        Optional evaluation grid vector.
    :type X: ndarray

    :returns:
        Psi - Weight matrix with Legendre functions up to ``degree``.
    :rtype: ndarray
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


def gb(degree, alpha, a=0, b=1, jmax=None, X=None):
    """Gegenbauer polynomials shifted to [a,b]

    :param degree:
        Degree of the polynomial.
    :type degree: int
    :param alpha:
        The Gegenbauer polynomials parameter.
    :type alpha: int
    :param a:
        Lower shift value (i.e., default - 0).
    :type a: int
    :param b:
        Upper shift value (i.e., default - 1).
    :type b: int
    :param jmax:
        Number of high-frequency lags.
    :param X:
        Optional evaluation grid vector.
    :type X: ndarray

    :returns:
        Psi - Weight matrix with Gegenbauer functions upto ``degree``.
    :rtype: ndarray
    """
    if jmax is not None:
        X = np.linspace(start=0, stop=1, num=jmax)

    if X is not None:
        raise ValueError("X is not provided. Either set X or set jmax.")

    n = len(X)
    P = np.ones([n, degree+2])
    Psi = np.ones([n, degree+1]) / math.sqrt(b-a)

    P[:, 1] = 2 * X / (b-a) - ((b+a)/(b-a))

    if degree > 0:
        for i in range(1, degree+1, 1):
            d = (2*i + 2*alpha + 1) * (2*i + 2*alpha + 2) / (2*(i+1)*(i + 2*alpha + 1))
            c = (alpha + 1)**2 * (2*i + 2*alpha + 2) / ((i+1) * (i + 2*alpha + 1) * (2*i + 2*alpha))
            P[:, i+1] = d * P[:, 1] * P[:, i] - c * P[:, i-1]
            Psi[:, i] = np.sqrt((2*i + 1) / (b-a) @ P[:, i])

    return Psi
