import numpy as np
from scipy.sparse.linalg import svds

def maxeig2(x):
    """
    Returns the largest squared singular value of a n-by-2 matrix x
    (the largest eigenvalue of the corresponding 2-by-2 matrix mat).
    """
    mat = np.dot(x.T, x)
    tr = mat[0, 0] + mat[1, 1]
    det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
    return (tr + np.sqrt(tr**2 - 4 * det)) / 2

def calc_gamma(x, ix, iy, bn):
    """
    Calculates a measure (gamma) for columns of matrix 'x' specified by ranges in 'ix' and 'iy'.
    """
    gamma = np.full(bn, np.nan)
    for g in range(bn):
        grabcols = slice(ix[g], iy[g] + 1)  # Python uses 0-based indexing
        submatrix = x[:, grabcols]
        ncols = submatrix.shape[1]
        
        if ncols > 2:
            # Calculate the largest singular value squared using RSpectra equivalent in scipy
            singular_values = svds(submatrix, k=1, return_singular_vectors=False)
            gamma[g] = singular_values[0]**2
        elif ncols == 2:
            gamma[g] = maxeig2(submatrix)
        else:
            gamma[g] = np.sum(submatrix**2)
    
    return gamma / x.shape[0]

