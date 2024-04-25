import numpy as np
import midasmlpy.src.sparsegl.sparsegllog_module as sgl # the sparse group lasso module from fortran
from scipy.sparse.linalg import svds
from sklearn.metrics import accuracy_score, roc_auc_score
from sklearn.model_selection import KFold

########################################################################################

########### Functions related to fitting the sparse group lasso model.##################

########################################################################################

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
            # Calculate the largest singular value squared 
            singular_values = svds(submatrix, k=1, return_singular_vectors=False)
            gamma[g] = singular_values[0]**2
        elif ncols == 2:
            gamma[g] = maxeig2(submatrix)
        else:
            gamma[g] = np.sum(submatrix**2)
    
    return gamma / x.shape[0]


def sgLasso_estimation(x, y, group_size, alsparse, pmax = 100, intr = True, nlam=None, ulam=None):
    """
    Implements the Sparse Group Lasso algorithm, which is a regularization technique combining 
    both lasso (L1) and group lasso (L2) penalties. This method is particularly useful for 
    scenarios where both individual feature sparsity and group sparsity are desired. The 
    algorithm assumes that features are divided into groups of equal size.

    Args:
        x (numpy.ndarray): A (n_obs, n_features) matrix of input features, where n_obs is the number of observations and 
                           n_features is the total number of features.
        y (numpy.ndarray): A vector of length n_obs that contains the dependent variable (response variable).
        group_size (int): The number of features in each group. It is assumed that all groups have the same number of features.
        alsparse (float): The alpha parameter that balances the L1 and L2 penalties. An alsparse close to 1.0 gives more weight to the 
                          L1 part (similar to Lasso), while an alsparse close to 0.0 gives more weight to the L2 part (similar to Group Lasso).
        nlam (int, optional): The number of lambda values to use for fitting the model. Default is 100.
        ulam (numpy.ndarray, optional): A sequence of lambda values to use for fitting the model. Default is np.ones(nlam).
        pmax (int, optional): The maximum number of non-zero coefficients allowed in the model. Default is 100.
        intr (bool, optional): If True, an intercept is included in the model. Default is True.

    Returns:
        tuple: Contains the following elements:
            nalam (int): The number of lambda values actually used.
            b0 (numpy.ndarray): The estimated intercept (if intr is True).
            beta (numpy.ndarray): Coefficient estimates for the predictors.
            activeGroup (numpy.ndarray): Indices of active groups in the model.
            nbeta (int): Number of active predictors.
            alam (numpy.ndarray): The sequence of lambda values used.
            npass (int): The number of passes (iterations) over the data.
            jerr (int): An error code (if any) from the fitting process (0 means no error).
    """
    if ulam is not None:
        ulam = np.asarray(ulam)
        nlam = len(ulam)  # Override nlam based on the length of ulam if ulam is provided
    elif nlam is None:
        nlam = 100  # Default value if neither ulam nor nlam is provided

    if ulam is None:
        ulam = np.ones(nlam)  # Default ulam if not provided

    nobs,nvars = x.shape[0], x.shape[1] # Number of observations and features
    eps = 1e-8 # Convergence threshold
    maxit = 3e8 # Maximum number of iterations
    bn = x.shape[1]//group_size # Number of groups as an integer
    bs = np.full(bn, group_size, dtype=int) # Elements in groups
    ix, iy =  list(range(0, nvars, group_size)), list(range(group_size-1, nvars, group_size)) # Placement og first column of each group in x
    gam = 0.25 * calc_gamma(x, ix, iy, bn) # Calculate gamma values for each group of features (columns) 
    pf, pfl1 = np.sqrt(bs),np.ones(nvars) # Penalty factors for L2 and L1 penalties
    dfmax = bn + 1 # Maximum number of groups
    ulam = np.ones(nlam) # Sequence of lambda values
    flmin = 0.01 if nobs < nvars else 1e-04
    lb,ub = np.full(bn, -np.inf),np.full(bn, np.inf) # Lower and upper bounds for the coefficients
    
    nalam, b0, beta, activeGroup, nbeta, alam, npass, jerr = sgl.log_sparse_four(x = x,
                    y = y, bn = bn, bs = bs, ix = ix, iy = iy, gam = gam, nobs = nobs, 
                    nvars = nvars, pf = pf, pfl1 = pfl1, dfmax = dfmax, pmax = pmax, 
                    nlam = nlam, flmin = flmin, ulam = ulam, eps = eps, maxit = maxit, 
                    intr = intr, lb = lb, ub = ub, alsparse = alsparse)
    return nalam, b0, beta, activeGroup, nbeta, alam, npass, jerr

########################################################################################

######### Functions related to finding the optimal sparse group lasso model.############

########################################################################################

def predict_binomial(x,b0,beta,threshold=0.5):
    """
    Predict binary outcomes (0 or 1) using logistic regression coefficients.

    Parameters:
    - x (ndarray): A 2D numpy array where each row represents a sample and each column represents a feature.
    - b0 (float): The intercept (bias) of the logistic regression model.
    - beta (ndarray): A 1D numpy array of coefficients for the logistic regression model corresponding to the features in `x`.
    - threshold (float, optional): The cutoff probability to determine the binary outcomes (defaults to 0.5).

    Returns:
    - predictions (ndarray): A 1D numpy array of binary outcomes (0 or 1). Each element corresponds to a sample in `x`.
    """
    probabilities = 1 / (1 + np.exp(-np.dot(x, beta) + b0))
    predictions = (probabilities > threshold).astype(int)
    return predictions

def evaluate_binomial(x, y, b0, beta,eval = 'auc', threshold=0.5):
    """
    Evaluate the performance of logistic regression models using specified metrics for different values of lambda.

    Parameters:
    - x (ndarray): A 2D numpy array of input features (identical to `x` in `predict`).
    - y (ndarray): A 1D numpy array containing the true binary outcomes (0 or 1) for each sample in `x`.
    - b0 (ndarray): A 1D numpy array of intercepts, one for each model being evaluated.
    - beta (ndarray): A 2D numpy array where each column corresponds to the coefficients of a model.
    - eval (str): The metric for evaluation; 'accuracy' for accuracy score, 'auc' for AUC score.

    Returns:
    - accuracies (list): If `eval` == 'accuracy', a list of accuracy scores for each model.
    - auc_scores (list): If `eval` == 'auc', a list of AUC scores for each model.
    """
    evaluation_score = [0] * len(b0)  # this will store evaluation score
    for l in range(len(b0)):
        predictions = predict_binomial(x,b0[l],beta[:,l], threshold=threshold)
        if eval == 'accuracy':
            evaluation_score[l] = accuracy_score(y, predictions)  
        elif eval == 'auc':
            evaluation_score[l] = roc_auc_score(y, predictions)
        else:
            raise ValueError("Invalid evaluation metric. Use 'accuracy' or 'auc'.")
    return evaluation_score

def cvlambda(x,y,alam,group_size, alsparse, pmax, intr,k_folds):
    """
    Perform k-fold cross-evaluation for logistic regression models 
    to evaluate the performance of each lambda.

    Parameters:
    - x (ndarray): A 2D numpy array of input features (identical to `x` in `predict`).
    - y (ndarray): A 1D numpy array containing the true binary outcomes (0 or 1) for each sample in `x`.
    - alam (ndarray): Array of lambda values to be cross validated over.
    - group_size (int): The number of features in each group.
    - alsparse (float): The alpha parameter that balances the L1 and L2 penalties.
    - pmax (int): The maximum number of non-zero coefficients allowed in the model.
    - k_folds (int, optional): The number of folds for cross-validation (defaults to 5).

    Returns:
    - mean_performace (ndarray): An array of mean AUC scores for each value of lambda.
    """
    # Split the data into k_folds
    kf = KFold(n_splits=k_folds)
    # initialize performance list
    performance = []
    for train_index, test_index in kf.split(x):
        # Based on the split, create the training and test data for this fold
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # Estimate the model on the training data
        _nalam, b0, beta, _activeGroup, _nbeta, _alam, _npass, _jerr = sgLasso_estimation(x_train, y_train, group_size, alsparse, pmax = pmax, intr = intr, ulam = alam)
        performance.append(evaluate_binomial(x_test, y_test, b0, beta,eval = 'auc', threshold=0.5))

    performance = np.array(performance)
    return np.mean(performance, axis=0) # return the mean performance across all folds

def bestmodel(x,y,group_size, alsparse, nlam = 100, pmax = 100, intr = True,k_folds = 5):
    """
    Find the best model using sparse group lasso. The sparse group lasso finds coefficients for nlam values of lambda, and the best model
    is chosen as the one with the highest mean performance in k-fold cross-validation.

    Parameters:
    - x (ndarray): A 2D numpy array of input features (identical to `x` in `predict`).
    - y (ndarray): A 1D numpy array containing the true binary outcomes (0 or 1) for each sample in `x`.
    - group_size (int): The number of lags used for the legendre polynomials.
    - alsparse (int): The balancing parameter alpha.
    - nlam (int, optional): The number of lambda values to evaluate (defaults to 100).
    - pmax (int, optional): The maximum number of non-zero coefficients allowed in the model (defaults to 100).
    - intr (bool, optional): If True, an intercept is included in the model (defaults to True).
    - k_folds (int, optional): The number of folds for cross-validation (defaults to 5).

    Returns:
    - best_model (dict): A dictionary containing the following
        - 'b0' (float): The intercept of the best model.
        - 'beta' (ndarray): The coefficients of the best model.
        - 'maximized_performance' (float): The maximized performance of the best model.
        - 'best_lambda' (int): The index of the best lambda value.
    """
    # Find best model
    _nalam, b0, beta, _activeGroup, _nbeta, alam, _npass, _jerr = sgLasso_estimation(x, y, group_size, alsparse, pmax, intr)

    # Find mean performance for each lambda
    mean_performance = cvlambda(x,y,alam,group_size, alsparse, pmax, intr,k_folds)
    best_lambda = np.argmax(mean_performance)
    max_performance = mean_performance[best_lambda]
    b0 = b0[best_lambda] 
    beta = beta[:,best_lambda] 
    return {
        "b0": b0,
        "beta": beta,
        "maximized_performance": max_performance,
        "best_lambda": best_lambda
    }
    