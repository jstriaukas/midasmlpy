# External imports
import numpy as np
from scipy.optimize import minimize


def optim_ardl_beta(y, z, x, poly_spec=0, nbtrials=100):
    """
    :param y:
    :param z:
    :param x:
    :param poly_spec:
    :param nbtrials:
    :return: coef
    :rtype: ndarray
    """
    n = len(y)
    k = None
    try:
        k = np.shape(z)[1]
    except IndexError:
        print("The input z for optim_ardl_beta has only one dimension, instead of the expected 2.")
    if k is not None:
        z = np.array(z).reshape((z, 1))
        k = 1

    # append a vector of ones
    iota = np.ones((n, 1))
    z = np.column_stack((iota, z))

    # store data into a dict
    data_dict = {"y": y, "z": z, "x": x, "poly_sec": poly_spec}
    # data_list = [y, z, x, poly_spec]

    # -------------------- main optimization --------------------#
    x0 = np.array([y, z, x, poly_spec, nbtrials])

    opt = minimize(estimate_ardl_beta, x0, method='Newton-CG',
                   jac=gradient_ardl_beta, hess=hessian_ardl_beta,
                   options={'xatol': 1.5e-10, 'disp': True})
    # -------------------- xxxxxxxxxxxxxxxx --------------------#

    # back-out parameters
    k1 = opt.args[0]
    k2 = opt.args[1]
    k3 = opt.args[2]
    beta = opt.args[3]
    c = opt.args[4]
    rho = opt.args[5][5:len(opt.args)]

    # sort the output
    other_x = None
    for i in range(1, k + 1):
        other_x = " ".join([other_x, "".join(["other_", i])])

    coef = None
    if poly_spec == 0:
        coef = np.array([*c, *rho, *beta, *k1, *k2, *k3]).reshape((1, k + 5))
    elif poly_spec == 1:
        coef = np.array([*c, *rho, *beta, *k2, *k3]).reshape((1, k + 4))
    elif poly_spec == 2:
        coef = np.array([*c, *rho, *beta, *k1, *k2]).reshape((1, k + 4))

    return coef


def get_constr_beta(poly_sec, k):
    pass


def estimate_ardl_beta(args, data_list):
    pass


def gradient_ardl_beta(args, data_list):
    pass


def hessian_ardl_beta(args, data_list):
    pass
