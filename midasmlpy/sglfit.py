import random
import numpy as np
import math
import matplotlib.pyplot as plt
from datetime import date
import sglfitF
import sys


Mypath = "C:/Users/Utilisateur/Documents/GitHub/midasmlpy/R_test"


class class_dgmatrix():
    def __init__(self):
        self.i = ""
        self.p = ""
        self.Dim = ""
        self.Dimnames = ""
        self.x = ""

class class_fit():
    def __init__(self):
        self.nalam = None
        self.b0 = None
        self.beta = None
        self.ibeta = None
        self.nbeta = None
        self.theta = None
        self.itheta = None
        self.ntheta = None
        self.alam = None
        self.npass = None
        self.jerr = None
        self.a0 = None
        self.nf = None
        
def sglfitpath(x, y, nlam, flmin, ulam, isd, intr, nf, eps, peps, dfmax, pmax, jd,
               pf, gindex, ngroups, maxit, gamma, nobs, nvars, vnames):
    # gamma setup
    if gamma < 0 or gamma > 1:
        sys.exit("gamma must be in [0,1]")
        
    # gindex should contain only the variable index at the end of each group - no overlap is allowed
    global ha
    ha = gindex
    gindex = find_diff(gindex)
    if any(i > 1 for i in gindex):
        sys.exit("only adjacent group memberships are allowed")
    gindex = find_index(gindex)
    ngroups = int(ngroups)
    #gindex = int(gindex)
    nobs = int(nobs)
    dfmax= int(dfmax)
    pmax = int(pmax)
    isd = int(isd)
    intr = int(intr)
    maxit = int(maxit)
    if nf == 0:   
        nalam, b0, beta, ibeta, nbeta, alam, npass, jerr = sglfitF.sglfit(
            gamma, ngroups, gindex, nobs, nvars, x, y, pf, dfmax,
            pmax, nlam, flmin, ulam, eps, peps, isd, 
            intr, maxit)
        nf = intr
        _ = np.transpose(beta)
        beta_ = []
        for i in range(len(_)):
            beta_ = beta_ + _[i].tolist()
        beta = np.array(beta_)

        fit = class_fit()        
        fit.nalam = nalam
        fit.b0 = b0
        fit.beta = beta
        fit.ibeta = ibeta
        fit.nbeta = nbeta
        fit.alam = alam
        fit.npass = npass
        fit.jerr = jerr
        fit.nf = nf
        
    else:
         print("not yet implemented")
    

    # output
    outlist = getoutput(fit, maxit, pmax, nvars, vnames)
    update_out = {"npasses": fit.npass, "jerr": fit.jerr}
    outlist.update(update_out)
    outlist["dimx"] = [nobs, nvars]
    #print("outlist:",outlist)

    return outlist
    
    
def getoutput(fit, maxit, pmax, nvars, vnames):
    nalam = fit.nalam
    nbeta = get_array_by_index(fit.nbeta, nalam)
    nbetamax = find_max(nbeta)
    lam = get_array_by_index(fit.alam, nalam)
    stepnames = get_step_names("s", nalam)
    n, errmsg = err(fit.jerr, maxit, pmax)
    if n == 1:
        sys.exit(errmsg)
    elif n == -1:
        print(errmsg)
    dd = [nvars, nalam]
    if nbetamax > 0:
        beta = find_beta(fit.beta, pmax, nalam, nbetamax)
        df_beta = sum(abs(beta) > 0)
        ja = get_array_by_index(fit.ibeta, nbetamax)
        oja = np.argsort(ja) + 1
        ja = ja[oja-1].tolist() * nalam
        ibeta = np.concatenate([[1], np.repeat(nbetamax, nalam)]).cumsum()
        
        _ = np.transpose(beta[(oja - 1).tolist()])
        __ = []
        for i in range(len(_)):
            __ = __ + _[i].tolist()
        
        beta_ = class_dgmatrix()
        beta_.i = MyOwnChange([k-1 for k in ja],"int")
        beta_.p = MyOwnChange([k-1 for k in ibeta],"int")
        beta_.x = __    
        beta_.Dimnames = [vnames,stepnames]
        beta_.Dim = dd
            
        beta = beta_
        
    else:
        beta = zeromat(nvars, nalam, vnames, stepnames)
        df_beta = [0] * nalam
    b0 = fit.b0
    
    if not b0 is None:
        b0 = get_array_by_index(b0, nalam)
    
    a0 = fit.a0
    nf = fit.nf
    
    if not a0 is None:
        a0 = np.full((nf, nalam), a0)
    
    if not fit.theta is None:
        ntheta = get_array_by_index(fit.ntheta, nalam)
        nthetamax = find_max(ntheta)
        if nthetamax > 0:
            theta = find_beta(fit.theta, pmax, nalam, nthetamax)
            df_theta = sum(abs(theta) > 0)
            ja = get_array_by_index(fit.itheta, nthetamax)
            oja = np.argsort(ja) + 1
            ja = ja[oja-1].tolist() * nalam
            itheta = np.concatenate([[1], np.repeat(nthetamax, nalam)]).cumsum()

            _ = np.transpose(theta[(oja - 1).tolist()])
            __ = []
            for i in range(len(_)):
                __ = __ + _[i].tolist()

            theta_ = class_dgmatrix()
            theta_.i = MyOwnChange([k-1 for k in ja],"int")
            theta_.p = MyOwnChange([k-1 for k in itheta],"int")
            theta_.x = __    
            theta_.Dimnames = [vnames,stepnames]
            theta_.Dim = dd
                
            theta = theta_
        else:
            theta = zeromat(nvars, nalam, vnames, stepnames)
            df_theta = [0] * nalam
            
        t0 = fit.t0
        if not t0 is None:
            t0 = get_array_by_index(t0, nalam)

        return {
                "b0": b0,
                "a0": a0,
                "beta": beta,
                "t0": t0,
                "theta": theta,
                "df_beta": df_beta,
                "df_theta": df_theta,
                "dd": dd,
                "lam": lam
                }
    return {
            "b0": b0,
            "a0": a0,
            "beta": beta,
            "df_beta": df_beta,
            "dd": dd,
            "lam": lam,
            "nf": nf
            }
    

def zeromat(nvars, nalam, vnames, stepnames):
    dgmatrix = class_dgmatrix()
    
    dgmatrix.i = [0]*nalam
    dgmatrix.p = list(range(nalam+1))
    dgmatrix.Dim = [nvars, nalam]
    dgmatrix.Dimnames = [vnames,stepnames]
    dgmatrix.x = [0]*nalam
    
    return dgmatrix

            
def find_beta(beta, pmax, nalam, nbetamax):
    beta = get_array_by_index(beta, pmax * nalam)
    beta = np.transpose(beta.reshape((nalam,pmax)))
    new_beta = get_array_by_index(beta, nbetamax)
    return new_beta

def get_array_by_index(x, index):
    arr = x[range(index)]
    return arr



def find_diff(x):
    diff_arr = []
    for i in range(len(x) - 1):
        diff_arr.append(x[i + 1] - x[i])
    return diff_arr



def find_index(x):
    index_arr = []
    for i in range(1,len(x)+1):
        if x[i-1] == 1:
            index_arr.append(i)
    index_arr.append(len(x)+1)
    return index_arr

def find_max(x):
    max_val = 0
    for i in range(len(x)):
        if x[i] > max_val:
            max_val = x[i]
    return max_val



def get_step_names(str1, step):
    name_arr = []
    for i in range(step):
        i = i
        names = '"' + str1 + str(i) + '"'
        name_arr.append(names)
    return name_arr

def err(n, maxit, pmax):
    if n == 0:
        msg = ""
    elif n > 0:
        if n < 7777:
            msg = "Memory allocation error"
        if n == 7777:
            msg = "All used predictors have zero variance"
        n = 1
        msg = "in fortran code -" + msg
    elif n < 0:
        if n > -10000:
            msg = "Convergence for " + str(-n) + "th lambda value not reached after maxit=" + str(
                maxit) + " iterations; solutions for larger lambdas returned "
        if n < -10000:
            msg = "Number of nonzero coefficients along the path exceeds pmax=" + pmax + " at" + str(
                -n - 10000) + "th lambda value; solutions for larger lambdas returned"
        n = -1
        msg = "from fortran code -" + msg
    
    
    return n, msg



with open(Mypath + '\\' + "input_x.txt") as datx:
    strx=datx.read()

xx=np.fromstring(strx, dtype=float, sep=' ')
x=np.reshape(xx,(100,20))

with open(Mypath + '\\' + "input_y.txt") as daty:
    stry=daty.read()
    
y=np.fromstring(stry, dtype=float, sep=' ')

def MyOwnChange(x, typee):
    if typee == "float":
        if (isinstance(x, list)):
            x = [float(z) for z in x]
        else:
            x = float(x)
    elif typee == "int":
        if (isinstance(x, list)):
            x = [int(z) for z in x]
        else:
            x = int(x)
    return x

def sglfit(x, y, gamma = 1.0, nlambda = 100, method = "single", nf = None,
           lambda_factor = None, lambda_ = None, pf = None, gindex = None,
           dfmax = None, pmax = None, standardize = False, 
           intercept = False, eps = 1e-08, maxit = 1000000, peps = 1e-08):
    
   
    np = x.shape
    nobs = np[0]
    nvars = np[1]
    
    if gindex == None:
        gindex = [i+1 for i in range(nvars)]
        
    ngroups = int(max(gindex))
    
    if lambda_factor == None:
        if nobs < nvars : 
            lambda_factor = 1e-02
        else:
            lambda_factor = 1e-04
    
    if pf == None:
        pf = [1 for i in range(nvars)]
    
    if dfmax == None:
        dfmax = nvars + 1
        
    if pmax == None:
        pmax = min(dfmax * 1.2, nvars)
        
    if maxit == None:
        maxit = 1000000
    
    if nvars > 0:
        vnames = ["V" + str(i+1) for i in range(nvars)]
    if len(y) != nobs:
        sys.exit("x and y have different number of observations") 
    
    if len(y.shape) > 1:
        sys.exit("Multivariate response is not supported now")
    
    if len(pf) != nvars:
        sys.exit("Size of L1 penalty factors does not match the number of input variables")
    if len(gindex) != nvars :
        sys.exit("Group index does not match the number of input variables")

    maxit = MyOwnChange(maxit,"int")
    pf = MyOwnChange(pf,"float")
    gindex = MyOwnChange(gindex,"int")
    isd = MyOwnChange(standardize,"int")
    intr = MyOwnChange(intercept,"int")
    eps = MyOwnChange(eps,"float")
    peps = MyOwnChange(peps,"float")
    dfmax = MyOwnChange(dfmax,"int")
    pmax = MyOwnChange(pmax,"int")
    jd = MyOwnChange(0,"int")
    
    if method == "single":
        nf = MyOwnChange(0,"int")
    else:
        isd = MyOwnChange(standardize,"int")
        if method == "pooled":
            intr = MyOwnChange(1,"int")
            nf = MyOwnChange(0,"int")
        elif method == "fe":
            intr = MyOwnChange(0,"int")
            if nf == None:
                sys.exit("Mehtod set as fixed effects without specifying the number of fixed effects (nf).")
            T = nobs/nf
            if round(T) == T:
                sys.exit("Number of fixed effects (nf) is not a multiple of the time series dimension, ie nf * T != nobs.")
    
    nlam = MyOwnChange(nlambda,"int")
    if lambda_ == None:
        if lambda_factor >= 1:
            sys.exit("lambda8factor should be less than 1")
        
        flmin = MyOwnChange(lambda_factor,"float")
        ulam = [float(0) for i in range(nlambda)]
        ulam[0] = -1
        ulam = MyOwnChange(ulam,"float")
    else:
        flmin = MyOwnChange(1,"float")
        if (isinstance(lambda_, list)):
            lambda_ = MyOwnChange(lambda_,"float")
            _ = [i for i in lambda_ if i < 0]
            if len(_) == 0:
                sys.exit("lambdas should be non-negative")
        else:
            lambda_ = MyOwnChange(lambda_,"float")
            if lambda_ < 0:
                sys.exit("lambda should be non-negative")
        
        lambda_.sort(reverse=True)
        
        ulam = MyOwnChange(lambda_,"float")
        nlam = MyOwnChange(len(lambda_),"int")
    
    fit = sglfitpath(x, y, nlam, flmin, ulam, isd, intr, nf, eps, peps, dfmax, pmax, jd, 
                     pf, gindex, ngroups, maxit, gamma, nobs, nvars, vnames)
    
    return fit
        


solution = sglfit(x, y, gamma = .5, nlambda = 100, method = "single",
                  gindex=[1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,4,4,4,4,4],
                  standardize = False, intercept = True)

