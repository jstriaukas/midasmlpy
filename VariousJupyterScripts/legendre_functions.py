import numpy as np

def legendre_polynomial(x, d):
    """
    This function calculates the Legendre polynomial of order n for the input x.
    For degrees 2 and above, the function uses the recursive relation to calculate the polynomial.
    """
    if d == 0:
        return 1
    elif d == 1:
        return x
    else:
        return ((2 * d - 1) * x * legendre_polynomial(x, d - 1) - (d - 1) * legendre_polynomial(x, d - 2)) / d

# Create matrix of legendre polynomials for lags and degrees
def legendre_matrix_create(x_lags, legendre_degree = 3):
    """
    This function creates a matrix of Legendre polynomials for the input x_lags and legendre_degree.
    """
    # Create list of equally spaced values between -1 and 1 for use in legendre polynomials
    x_values = np.linspace(-1, 1, num=x_lags)
    # Loop through the lags and degrees to create a lag x degree matrix of legendre polynomials
    legendre_matrix = [[legendre_polynomial(x, n) for n in range(legendre_degree)] for x in x_values]
    return np.array(legendre_matrix)