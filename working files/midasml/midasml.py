import midasml2.midasml2 as ml


def lin_grad_calc(_nrow, _eta, _Y, _ldot):
    ml.lin_grad_calc(_nrow, _eta, _Y, _ldot)
    return _ldot


def lin_neg_log_likelihood_calc(_nrow, _eta, _y):
    return ml.lin_neg_log_likelihood_calc(_nrow, _eta, _y)


def lin_solver_lam_grid(
        _X, _Y, _index, _nrow, _ncol, _numGroup, _beta,
        _rangeGroupInd, _groupLen, _lambda1, _lambda2,
        _innerIter, _thresh, _ldot, _nullBeta, _gamma, _eta,
        _betaIsZero, _groupChange, _isActive, _useGroup,
        _step, _reset
):
    ml.lin_solver_lam_grid(
        _X, _Y, _index, _nrow, _ncol, _numGroup,
        _beta, _rangeGroupInd, _groupLen, _lambda1,
        _lambda2, _innerIter, _thresh, _ldot, _nullBeta,
        _gamma, _eta, _betaIsZero, _groupChange, _isActive, _useGroup,
        _step, _reset
    )


def cvfolds(_nfolds, _nrow):
    return ml.cvfold_s(_nfolds, _nrow)


def getmin(_lambda, _cvm, _cvsd, _which_lambda):
    return ml.getmin(_lambda, _cvm, _cvsd, _which_lambda)


def sgl_fit(_Z, _Y, _lambda1, _lambda2, _innerIter, _outerIter, _thresh, _outerThresh, _gamma_solver, _step, _reset):
    return ml.sgl_fit(_Z, _Y, _lambda1, _lambda2, _innerIter, _outerIter, _thresh, _outerThresh, _gamma_solver, _step,
                      _reset)


def sgl_fit_path(_X, _Z,
                 _y, _index, _dummies,
                 _l1_frac, _l21_frac,
                 _dummies_index, _lambdas,
                 _gamma_w, _innerIter, _outerIter, _thresh, _outerThresh, _gamma_solver, _step, _reset):
    return ml.sgl_fit_path(
        _X, _Z, _y, _index, _dummies, _l1_frac, _l21_frac,
        _dummies_index, _lambdas, _gamma_w, _innerIter, _outerIter, _thresh, _outerThresh,
        _gamma_solver, _step, _reset
    )


def fastols(_Y, _X, _intercept):
    return ml.fastol_s(_Y, _X, _intercept)


def fastals(_Y, _X, _tau, _intercept, _maxIter, _thresh):
    return ml.fastal_s(_Y, _X, _tau, _intercept, _maxIter, _thresh)


def fastrq(_Y, _X, _intercept, _tau):
    return ml.fastr_q(_Y, _X, _intercept, _tau)


def midas_pr(_Y, _X, _intercept, _tau, _which_loss, _num_evals):
    return ml.midaspr(_Y, _X, _intercept, _tau, _which_loss, _num_evals)


def midasar_pr(_Y, _YLAG, _X, _intercept, _tau, _which_loss, _num_evals):
    return ml.midasarpr(_Y, _YLAG, _X, _intercept, _tau, _which_loss, _num_evals)
