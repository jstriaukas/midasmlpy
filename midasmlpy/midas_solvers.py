# External imports
from multiprocessing import Pool
import numpy as np
import scipy


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
    x0 = np.array([y, z, x, poly_spec, nbtrials])  # with same info as contained in data_dict

    def multistart(num_core=2, count=100):
        async_res = []
        with Pool(processes=num_core) as p:
            for _ in range(count):
                async_res.append(p.apply_async(
                    scipy.optimize.minimize(estimate_ardl_beta, x0, method='Newton-CG',
                                            jac=gradient_ardl_beta, hess=hessian_ardl_beta,
                                            options={'xatol': 1.5e-10, 'disp': True})
                ))
        return min(async_res)

    opt = multistart()

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


def estimate_ardl_beta(args, data_dict):
    pass


def gradient_ardl_beta(args, data_dict):
    x = data_dict['x']
    y = data_dict['y']
    z = data_dict['z']

    p = None
    try:
        p = np.shape(x)[1]
    except IndexError:
        print("The input x for gradient_ardl_beta has only one dimension, instead of the expected 2.")

    ii = np.ones([p, 1])
    xx = np.arange(1, p + 1) / (p + 1)

    poly_spec = data_dict['poly_spec']

    if poly_spec == 0:
        theta1 = args[0]
        theta3 = args[2]
    elif poly_spec == 1:
        theta1 = 1
        theta3 = args[2]
    elif poly_spec == 2:
        theta1 = args[0]
        theta3 = 0
    elif poly_spec == 3:
        theta1 = 1
        theta3 = 0

    theta2 = args[1]
    beta = args[3]
    rho = args[4:]
    m = scipy.special.gamma(theta1+theta2/(scipy.special.gamma(theta1) * scipy.special.gamma(theta2)))
    weights = xx**(theta1 - 1) * (ii - xx)**(theta2 - 1) + theta3

    nabla_G = 2 * (y - beta * np.matmul(x, weights) - np.matmul(z, rho))

    T_star_w_C = -beta * x.T
    T_star_beta_C = np.matmul(-1 * weights.T, x.T)
    T_star_rho_C = -1 * z.T

    weights_tilde = xx**(theta1 - 1) * (ii - xx)**(theta2 - 1) * m

    log_derivative_theta2 = np.log(ii - xx) + scipy.special.digamma(theta1 + theta2) - scipy.special.digamma(theta2)
    T_star_theta2_W = np.transpose(log_derivative_theta2*weights_tilde)

    part2 = np.matmul(np.matmul(T_star_theta2_W, T_star_w_C), nabla_G)  # theta2

    if poly_spec == 0:
        log_derivative_theta1 = np.log(xx) + scipy.special.digamma(theta1 + theta2) - scipy.special.digamma(theta1)
        T_star_theta1_W = np.transpose(log_derivative_theta1 * weights_tilde)
        T_star_theta3_W = ii.T
        part1 = np.matmul(np.matmul(T_star_theta1_W, T_star_w_C), nabla_G)  # theta1
        part3 = np.matmul(np.matmul(T_star_theta3_W, T_star_w_C), nabla_G)  # theta3
    elif poly_spec == 1:
        T_star_theta3_W = ii.T
        part1 = 0
        part3 = np.matmul(np.matmul(T_star_theta3_W, T_star_w_C), nabla_G)  # theta3
    elif poly_spec == 2:
        log_derivative_theta1 = np.log(xx) + scipy.special.digamma(theta1 + theta2) - scipy.special.digamma(theta1)
        T_star_theta1_W = np.transpose(log_derivative_theta1 * weights_tilde)
        part1 = np.matmul(np.matmul(T_star_theta1_W, T_star_w_C), nabla_G)  # theta1
        part3 = 0
    elif poly_spec == 3:
        part1 = 0
        part3 = 0

    part4 = np.matmul(T_star_beta_C, nabla_G)  # beta
    part5 = np.matmul(T_star_rho_C, nabla_G)  # rho
    grad = np.concatenate([part1, part2, part3, part4, part5]).flatten()

    return grad


# TODO: Complete the following
def hessian_ardl_beta(args, data_dict):
    pass
