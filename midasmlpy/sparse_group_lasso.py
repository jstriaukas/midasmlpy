import numpy as np
from scipy.sparse.linalg import svds
from sklearn.metrics import accuracy_score, roc_auc_score, mean_squared_error, r2_score
from sklearn.model_selection import StratifiedKFold, KFold
import random

from midasmlpy.src.sparseglf90 import sparsegllog_module

random.seed(111)


########################################################################################

########### Functions related to fitting the sparse group LASSO model.##################

########################################################################################

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
            singular_values = svds(submatrix, k=1, return_singular_vectors=False, random_state=42)
            gamma[g] = singular_values[0] ** 2
        elif ncols == 2:
            # Returns the largest squared singular value of a n-by-2 matrix x
            mat = np.dot(submatrix.T, submatrix)
            tr = mat[0, 0] + mat[1, 1]
            det = mat[0, 0] * mat[1, 1] - mat[0, 1] * mat[1, 0]
            gamma[g] = (tr + np.sqrt(tr ** 2 - 4 * det)) / 2
        else:
            gamma[g] = np.sum(submatrix ** 2)

    return gamma / x.shape[0]


def sgLASSO_estimation(x, y, group_size, alsparse, family='binomial', pmax=100, intr=True, nlam=None, ulam=None):
    """
    Implements the sparse group LASSO algorithm, which is a regularization technique combining 
    both LASSO (L1) and group LASSO (L2) penalties. This method is particularly useful for 
    scenarios where both individual feature sparsity and group sparsity are desired. The 
    algorithm assumes that features are divided into groups of equal size.

    Args:
        x (numpy.ndarray): A (n_obs, n_features) matrix of input features, where n_obs is the number of observations and 
                           n_features is the total number of features.
        y (numpy.ndarray): A vector of length n_obs that contains the dependent variable (response variable).
        group_size (int): The number of features in each group. It is assumed that all groups have the same number of features.
        alsparse (float): The alpha parameter that balances the L1 and L2 penalties. An alsparse close to 1.0 gives more weight to the 
                          L1 part (similar to LASSO), while an alsparse close to 0.0 gives more weight to the L2 part (similar to Group LASSO).
        family (str): The type of model to fit. Possible values are 'binomial' for logistic regression and 'gaussian' for linear regression.
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
        ulam = np.array(ulam)
        nlam = len(ulam)  # Override nlam based on the length of ulam if ulam is provided
    elif nlam is None:
        nlam = 100  # Default value if neither ulam nor nlam is provided

    if ulam is None:
        ulam = np.ones(nlam)  # Default ulam if not provided

    nobs, nvars = x.shape[0], x.shape[1]  # Number of observations and features
    eps = 1e-8  # Convergence threshold
    maxit = 1000000  # Maximum number of iterations
    bn = x.shape[1] // group_size  # Number of groups as an integer
    bs = np.full(bn, group_size, dtype=int)  # Elements in groups
    ix, iy = np.array(range(0, nvars, group_size)), np.array(
        range(group_size - 1, nvars, group_size))  # Placement og first column of each group in x
    pf, pfl1 = np.sqrt(bs), np.ones(nvars)  # Penalty factors for L2 and L1 penalties
    dfmax = bn + 1  # Maximum number of groups
    flmin = 0.01 if nobs < nvars else 1e-04
    lb, ub = np.full(bn, -np.inf), np.full(bn, np.inf)  # Lower and upper bounds for the coefficients

    jerr = 0
    npass = 0
    b0 = None
    beta = np.zeros(nvars)
    alam = np.zeros(nlam)
    mse = None

    if family == 'binomial':
        gam = 0.25 * calc_gamma(x, ix, iy, bn)  # Calculate gamma values for each group of features (columns)
        _nalam, b0, beta, _activeGroup, _nbeta, alam, npass, jerr = (sparsegllog_module.log_sparse_four(x=x,
                                                                                              y=y, bn=bn,
                                                                                              bs=bs,
                                                                                              ix=ix + 1,
                                                                                              iy=iy + 1,
                                                                                              # iy and ix are +1 as fortran is index 1 while python is index 0
                                                                                              gam=gam,
                                                                                              nobs=nobs,
                                                                                              nvars=nvars,
                                                                                              pf=pf,
                                                                                              pfl1=pfl1,
                                                                                              dfmax=dfmax,
                                                                                              pmax=pmax,
                                                                                              nlam=nlam,
                                                                                              flmin=flmin,
                                                                                              ulam=ulam,
                                                                                              eps=eps,
                                                                                              maxit=maxit,
                                                                                              intr=intr,
                                                                                              lb=lb, ub=ub,
                                                                                              alsparse=alsparse))
        mse = None  # to make it easier to return the same number of variables for all families
    if family == 'gaussian':
        if intr:
            y = y - y.mean()
        gam = calc_gamma(x, ix, iy, bn)  # Calculate gamma values for each group of features (columns)
        _nalam, b0, beta, _activeGroup, _nbeta, alam, npass, jerr, mse = sparsegllog_module.sparse_four(x=x,
                                                                                                     y=y, bn=bn,
                                                                                                     bs=bs,
                                                                                                     ix=ix + 1,
                                                                                                     iy=iy + 1,
                                                                                                     # iy and ix are +1 as fortran is index 1 while python is index 0
                                                                                                     gam=gam,
                                                                                                     nobs=nobs,
                                                                                                     nvars=nvars,
                                                                                                     pf=pf,
                                                                                                     pfl1=pfl1,
                                                                                                     dfmax=dfmax,
                                                                                                     pmax=pmax,
                                                                                                     nlam=nlam,
                                                                                                     flmin=flmin,
                                                                                                     ulam=ulam,
                                                                                                     eps=eps,
                                                                                                     maxit=maxit,
                                                                                                     intr=intr,
                                                                                                     lb=lb, ub=ub,
                                                                                                     alsparse=alsparse)
    if jerr != 0:
        raise ValueError("Error in the sparse group LASSO estimation.")
    if npass == maxit:
        raise ValueError("Failed to converge in the sparse group LASSO estimation.")
    return b0, beta, alam, npass, jerr, mse


########################################################################################

######### Functions related to finding the optimal sparse group LASSO model. ###########

########################################################################################

def predict_binomial(x, b0, beta, threshold=0.5):
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
    if b0 is None:
        b0 = 0

    probabilities = 1 / (1 + np.exp(-np.dot(x, beta) + b0))
    predictions = (probabilities > threshold).astype(int)
    return predictions


def evaluate_binomials(x, y, b0, beta, eval='auc', threshold=0.5):
    """
    Evaluate the performance of several logistic regression models using specified metrics for different values of lambda.

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
        predictions = predict_binomial(x, b0[l], beta[:, l], threshold=threshold)
        if eval == 'accuracy':
            evaluation_score[l] = accuracy_score(y, predictions)
        elif eval == 'auc':
            evaluation_score[l] = roc_auc_score(y, predictions)
        else:
            raise ValueError("Invalid evaluation metric. Use 'accuracy' or 'auc'.")
    return evaluation_score


def predict_gaussian(x, b0, beta):
    """
    Predict gaussian outcomes using linear regression coefficients.

    Parameters:
    - x (ndarray): A 2D numpy array where each row represents a sample and each column represents a feature.
    - b0 (float): The intercept (bias) of the logistic regression model.
    - beta (ndarray): A 1D numpy array of coefficients for the logistic regression model corresponding to the features in `x`.

    Returns:
    - predictions (ndarray): A 1D numpy array of outcomes. Each element corresponds to a sample in `x`.
    """
    predictions = np.dot(x, beta) + b0
    return predictions


def evaluate_gaussian(x, y, b0, beta, intr, eval='mse'):
    """
    Evaluate the performance of several linear regression models using specified metrics.

    Parameters:
    - x (ndarray): A 2D numpy array of input features (identical to `x` in `predict`).
    - y (ndarray): A 1D numpy array containing the true continuous outcomes for each sample in `x`.
    - b0 (ndarray): A 1D numpy array of intercepts, one for each model being evaluated.
    - beta (ndarray): A 2D numpy array where each column corresponds to the coefficients of a model.
    - eval (str): The metric for evaluation; options are 'mse' for Mean Squared Error or 'r2' for R-squared.

    Returns:
    - evaluation_scores (list): A list of scores, either MSE or R-squared, for each model.
    """
    evaluation_scores = [0] * len(b0)  # this will store evaluation scores
    for l in range(len(b0)):

        predictions = predict_gaussian(x, b0[l], beta[:, l])
        if intr:  # Adjust predictions if intercept was used during fitting
            predictions += y.mean()
        if eval == 'mse':
            evaluation_scores[l] = mean_squared_error(y, predictions)
        elif eval == 'r2':
            evaluation_scores[l] = r2_score(y, predictions)
        else:
            raise ValueError("Invalid evaluation metric. Use 'mse' or 'r2'.")
    return evaluation_scores


def best_lambda_find(x, y, group_size, alsparse, family='binomial', nlam=100, pmax=100, intr=True, k_folds=5):
    """
    Find the best model using sparse group LASSO. The sparse group LASSO finds coefficients for nlam values of lambda, and the best model
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
        - 'best_lambda' (float): The lambda value of the best model.
    """
    # Find model nlam number of models
    b0, beta, alam, _npass, _jerr, mse = sgLASSO_estimation(x, y, group_size, alsparse, family, pmax, intr)

    kf = StratifiedKFold(n_splits=k_folds, shuffle=False, random_state=42)

    # Find mean performance for each lambda
    # Split the data into k_folds
    if family == 'binomial':
        kf = StratifiedKFold(n_splits=k_folds, shuffle=False, random_state=42)
    if family == 'gaussian':
        kf = KFold(n_splits=k_folds, shuffle=False, random_state=42)

    best_lambda = None

    # initialize performance list
    performance = []
    for train_index, test_index in kf.split(x, y):
        # Based on the split, create the training and test data for this fold
        x_train, x_test = x[train_index], x[test_index]
        y_train, y_test = y[train_index], y[test_index]
        # Estimate the model on the training data
        b0test, betatest, _alam, _npass, _jerr, msetrain = sgLASSO_estimation(x_train, y_train, group_size, alsparse,
                                                                              family, pmax=pmax, intr=intr, ulam=alam)
        if family == 'gaussian':
            performance.append(evaluate_gaussian(x_test, y_test, b0test, betatest, intr, eval='mse'))
        if family == 'binomial':
            performance.append(evaluate_binomials(x_test, y_test, b0test, betatest, eval='auc', threshold=0.5))

    performance = np.array(performance)
    mean_performance = np.mean(performance, axis=0)
    if family == 'gaussian':
        best_lambda = np.argmin(mean_performance)
    if family == 'binomial':
        best_lambda = np.argmax(mean_performance)
    return {'b0': b0,
            'beta': beta[:, best_lambda],
            'best_performance': mean_performance[best_lambda],
            'best_lambda': alam[best_lambda]}


def best_model(x, y, group_size, family='binomial', nlam=100, pmax=100, intr=True, k_folds=5, disp_flag=True, alpha_values=None,
               alpha=None):
    """
    Function to find the best model based on the maximized performance of the model. The function uses the bestlambda function to find the best lambda value for the model.

    Parameters:
    x: numpy array
        The predictors for the model
    y: numpy array
        The target variable for the model
    group_size: int
        The number of groups in the model
    nlam: int
        The number of lambda values to test
    pmax: int
        The maximum number of variables in the model
    intr: bool
        Whether to include the intercept in the model
    k_folds: int
        The number of folds to use in the cross-validation
    disp_flag: bool
        Whether to display the performance at different values of alpha
    alpha_values: int
        The number of alpha values to test
    alpha: list
        The alpha value to test

    Returns:
    dict
        A dictionary containing the best alpha value, the best performance, the intercept and the coefficients of the model
    """

    best_performance = None
    best_alsparse = None
    b0, beta = None, None  # Initialize parameters that will store best model coefficients

    if alpha is not None:
        alsparse_values = np.array(alpha)
    elif alpha_values is not None:
        alsparse_values = np.linspace(1, 0, alpha_values)
    else:
        alsparse_values = np.linspace(1, 0, 5)

    # Dictionary to store the average maximized performances for each alsparse
    if disp_flag:
        performance_dict = {}

    # Cross-validation process
    for alsparse in alsparse_values:
        model_result = best_lambda_find(x, y, group_size, alsparse, family, nlam=nlam, pmax=pmax, intr=intr,
                                        k_folds=k_folds)
        # Append the maximized performance of this fold
        if disp_flag:
            performance_dict[alsparse] = model_result['best_performance'].round(5)

        # If this fold has a higher maximized performance than the previous best, update the best performance
        if best_performance is None:
            best_performance = model_result['best_performance']
            best_alsparse = alsparse
            b0 = model_result['b0']
            beta = model_result['beta']
            best_lambda = model_result['best_lambda']
        else:
            if family == 'gaussian':
                if model_result['best_performance'] < best_performance:
                    best_performance = model_result['best_performance']
                    best_alsparse = alsparse
                    b0 = model_result['b0']
                    beta = model_result['beta']
                    best_lambda = model_result['best_lambda']
            if family == 'binomial':
                if model_result['best_performance'] > best_performance:
                    best_performance = model_result['best_performance']
                    best_alsparse = alsparse
                    b0 = model_result['b0']
                    beta = model_result['beta']
                    best_lambda = model_result['best_lambda']

    if disp_flag:
        print('The performance at different values of alpha are:')
        print(performance_dict)

    return {'best_alsparse': best_alsparse,
            'best_performance': best_performance,
            'b0': b0,
            'beta': beta,
            'best_lambda': best_lambda}
