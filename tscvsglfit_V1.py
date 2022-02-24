import random
import numpy as np
import math
from datetime import date


def tscv_sglfit(x1, y1, nlambda, method, nf, lamb, gindex, l, parallel, seed, standardize, intercept):
    N = x1.shape[0]
    p = x1.shape[1]
    y = y1
    sglfit_object = sglfit(x1, y1, nlambda, method, nf, lamb, gindex, standardize, intercept)

    # ---------------------------------------------------------
    lamb = sglfit_object['lamb']
    # predict -> coef
    # nz =
    if l <= 1:
        print("l must be at least 2; l=5 recommended")
    if K < 1:
        print("K must be at least 1; K=20 recommended")
    if K >= N:
        print("K must be at most T-1; K=20 recommended")
    if seed is None:
        seed = (date.today() - date(1970, 1, 1)).days
    random.seed(seed)
    foldid = [random.choice(range(N)) for _ in range(K)]
    outlist = np.empty((K, 1))
    if parallel:
        for i in range(K):
            whichfoldnot = foldid[i]
            whichgaptrain = computegapobs(whichfoldnot, N, l)
            # -----------------------------------------------
            # find whichgaptrain
            y_sub = find_whichgaptrain(y, whichgaptrain)
            x_sub = find_whichgaptrain(x, whichgaptrain)
            outlist[i] = sglfit(x_sub, y_sub, lamb, gamma, gindex, 'single')
    else:
        for i in range(K):
            whichfoldnot = foldid[i]
            whichgaptrain = computegapobs(whichfoldnot, N, l)
            # -----------------------------------------------
            # find whichgaptrain
            y_sub = find_whichgaptrain(y, whichgaptrain)
            x_sub = find_whichgaptrain(x, whichgaptrain)
            outlist[i] = sglfit(x_sub, y_sub, lamb, gamma, gindex, 'single')
    cvstuff = tscv_sglpath(outlist, lamb, x, y, foldid)

    tscv_sglpath(outlist, lamb, x, y, foldid)


def tscv_sglpath(outlist, lamb, x, y, foldid):
    typenames = 'Single outcome sg-LASSO'
    y = y
    K = len(foldid)
    # predmat = np.empty((len(y), len(lamb)))
    nlams = [0] * K
    for i in range(K):
        whichfold = foldid[i]
        fitobj = outlist[i]
        # ------------------------------------------
        # preds =
        # ------------------------------------------
        nlami = len(outlist[i]['lamb'])
        # predmat(whichfold, ) = preds
        nlams[i] = nlami
    cvraw = ()
    print(nlams)


def sglfit(x2, y2, nlambda, method, nf, lamb, gindex, standardize, intercept):
    # data setup
    method = method
    np = x2.shape
    nobs = int(np[0])
    nvars = int(np[1])
    vnames = ["V1", "V2", "V3", "V4", "V5", "V6", "V7", "V8", "V9", "V10", "V11", "V12", "V13", "V14", "V15", "V16",
              "V17", "V18", "V19", "V20"]
    ngroups = int(find_max(gindex))

    # parameter setup
    maxit = 1000000
    pf = [1] * nvars
    isd = int(standardize)
    intr = int(intercept)
    eps = 0.00000001
    peps = 0.00000001
    dfmax = int(nvars + 1)
    pmax = int(0)

    # panel regression setup
    if method == 'single':
        nf = int(0)
    else:
        isd = int(standardize)
        if method == 'pooled':
            intr = int(1)
            nf = int(0)
        elif method == 'fe':
            intr = int(0)
            if nf is None:
                print("Mehtod set as fixed effects without specifying the number of fixed effects (nf).")
            T = nobs / nf
            if round(T) is not T:
                print("Number of fixed effects (nf) is not a multiple of the time series dimension, ie nf * T != nobs.")

    # lambda setup
    nlam = int(nlambda)
    if nobs < nvars:
        lamb_factor = 0.01
    else:
        lamb_factor = 0.0001
    if lamb is None:
        if lamb_factor >= 1:
            print("lambda.factor should be less than 1")
        flmin = lamb_factor
        ulam = [0] * nlambda
        ulam[0] = -1
    else:
        flmin = [0] * 1
        if any(i >= 0 for i in lamb):
            print("lambdas should be non-negative")
        ulam = [0] * nlambda
        ulam[0] = -1
        nlam = ulam
        ulam = ulam.sort()

    #################################################################################
    fit = sglfitpth(x, y, ulam, intr, nf, gindex, gamma)


def sglfitpth(x, y, ulam, intr, nf, gindex, gamma):
    # gamma setup
    if gamma < 0 or gamma > 1:
        print("gamma must be in [0,1]")
    # gindex should contain only the variable index at the end of each group - no overlap is allowed
    gindex = find_diff(gindex)
    if any(i > 1 for i in gindex):
        print("only adjacent group memberships are allowed")
    gindex = find_index(gindex)

    # call Fortran
    if nf == 0:
        fit = fortran_sglfitF(gamma, ulam)

    # output
    if nf == 0:
        nf = intr


def fortran_sglfitF(gamma, ngroups, gindex, nobs, nvars, x, y, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, isd, intr, maxit, nalam, b0, beta, ibeta, nbeta, alam, npass, jerr):
    maxlam = 0
    if ulam[0] == -1:
        maxlambda(nvars, nobs, x, y, gamma, gindex, ngroups, pf, maxlam)
        for j in range(1, len(ulam)):
            tmp = math.log(maxlam) + (math.log(maxlam*flmin)-math.log(maxlam)*(j - 1)/(nlam - 1))
            ulam[j] = math.exp(tmp)
        ulam[0] = maxlam
    if gamma == 1:
        fortran_lassofitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, ibeta, nbeta, alam, npass, jerr, intr)
    else:
        fortran_sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, ibeta, nbeta, alam, npass, jerr, intr)
    if jerr > 0:
        return
    for l in range(len(nalam)):
        nk = nbeta(0)
        if isd == 1:
            for j in range(len(nk)):
                beta(j, l) = beta(j, l)/ibeta(j)
        if inr == 1:
            b0[0] = b0[0] - np.multiply(beta)

def fortran_sglfitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, ibeta, nbeta, alam, npass, jerr, intr):
    b = 0
    oldbeta = 0
    m = 0
    mm = 0
    npass = 0
    ni = npass
    mnl = min(mnlam, nlam)
    ju = 0
    r = y
    for k in range(len(ngroups)):
        gend = gindex[k]
        if k == 1:
            gstart = 1
        else:
            gstart = gindex[k-1] + 1
        tmp = sum()
        steps[k] = 1/math.sqrt(tmp)
    #------------------lambda loop(outmost loop)----------------#
    for l in range(len(nlam)):
        al = ulam[l]
        ctr = 0
        pln = 0
        #---------------outer loop--------------------------#
        npass = npass + 1
        dif = 0
        if intr == 1:
            oldbeta[0] = b[0]
        if ni >0:
            oldbeta[m[1:ni]] = b[m[1:ni]]
        pln = pln + 1
        for k in range(len(ngroups)):
            gend = gindex(k)
            if k == 1:
                gstart = 1
            else:
                gstart = gindex[k-1] + 1
            gs = gend - gstart+1
            gw = 0
            for gj in range(gstart, gend):
                gw = gw + pf(gj)
            gw = math.sqrt(gw)
            skip = 1
            if pln == 1:
                skip = 0
            if ju[k] == 1:
                skip = 0
            if skip == 0:
                #--------------sg-lasso proximal map-------------------#
                fortran_prox_sgl(gstart, gend, nvars, nobs, x, r, b[1:nvars], al, gamma, pf, peps, gw, steps[k])
                #--------------update remaining variables---------------#
                for g in range(gstart, gend):
                    if abs(b[g]>0):
                        if pln == 1:
                            ju[k] = 1
                        d = oldbeta[g]-b[g]
                        dif = max(dif, maj[g]*d^2)
                        if mm[g] == 0:
                            ni = ni + 1
                            if ni> pmax:
                                break
                            mm[g] = ni
                            m[ni] = g
        if ni >pmax:
            break
        if intr == 1:
            oldb = oldbeta[0]
            u = 0
            d = sum(r)/nobs
            if d^2 < eps:
                exit()
            b[0] = b[0] +d
            r = r-d
            d = b[0]-oldb
            if abs(d) > 0:
                dif = max(dif, d^2)
        if ni > pmax:
            exit()
        vrg = 1
        if intr == 1:
            if (b[0]-oldbeta[0])^2 >= eps:
                vrg = 0
        for j in range(1, ni):
            if(b[m[j]]-oldbeta[m[j]])^2 >= eps:
                vrg = 0
        if vrg == 1:
            exit()
        ctr = ctr + 1
        if ctr >maxit :
            jerr = -l
    #------------------final update & save results-------------#
    if ni > pmax:
        jerr = -10000 - l
        exit()
    if ni>0:
        beta[1:ni] = b[m[1:ni]]
        nbeta[l] = ni
        b0[l] = al
        alam[l] = l
        if flmin > 1:
            me = find_count(abs(beta[1:ni])> 0)
        if me > dfmax:
            exit()

def fortran_lassofitpathF(maj, gamma, ngroups, gindex, nobs, nvars, x, y, ju, pf, dfmax, pmax, nlam, flmin, ulam, eps, peps, maxit, nalam, b0, beta, ibeta, nbeta, alam, npass, jerr, intr):
    #----------------initialization------------------#
    b = 0.0
    oldbeta = 0.0
    m = 0
    mm = 0
    npass = 0
    ni = npass
    mnl = min(mnlam, nlam)
    ju = 0
    r = y


    # ------------------lambda loop(outmost loop)----------------#
    for l in range(len(nlam)):
        al = ulam[l]
        ctr = 0
        # ---------------outer loop--------------------------#
        if intr == 1:
            oldbeta[0] = b[0]
        if ni>0:
            oldbeta[m[1:ni]] = b[m[1:ni]]
        pln = 0
        npass = npass + 1
        dif = 0.0
        for k in range(1, nvars):
            if pln == 1:
                oldb = b[k]
                u = maj[k]*b[k] + np.multiply(r, x[:, k]/nobs)
                v = abs(u) - al*pf[k]
                if v>0.0:
                    tmp = math.copysign(v, u)/maj[k]
                else:
                    tmp = 0.0
                d = tmp - b[k]
                if d^2 <eps:
                    exit()
                b[k] = tmp
            d = b[k] - oldb
            if abs(d) > 0.0:
                dif = dif + maj[k] *d^2
                if mm[k] == 0:
                    ni = ni + 1
                    if ni> pmax:
                        exit()
                    mm[k] = ni
                    m[ni] = k
            if abs(b[k])>0.0:
                ju[k] = 1
        else:
            if ju[k] == 1:
                oldb = b[k]
                u = np.multiply(r, x[:, k])
                u = maj[k] * b[k] + u/nobs
                v = abs(u) - al*pf[k]
                if v>0.0:
                    tmp = math.copysign(v, u)/maj[k]
                else:
                    tmp = 0.0
                d = tmp - b[k]
                if d^2 < eps:
                    exit()
                b[k] = tmp
                r = r - x[:, k]*d
            d = b[k]-oldb
            if abs(d)>0.0:
                dif = max(dif, maj[k]*d^2)
        if ni>pmax:
            exit()
        if intr == 1:
            oldb = b[0]
            d = sum(r)/nobs
            if d^2 < eps:
                exit()
            b[0] = b[0] + d
            r = r -d
            d = b[0] - oldb
            if abs(d)> 0.0:
                dif = max(dif, d^2)
            if dif <eps:
                exit()
        if ni >pmax:
            exit()
        vrg = 1

        if intr == 1:
            if (b[0] - oldbeta[0]) ^ 2 >= eps:
                vrg = 0
        for j in range(1, ni):
            if (b[m[j]] - oldbeta[m[j]]) ^ 2 >= eps:
                vrg = 0
        if vrg == 1:
            exit()
        ctr = ctr + 1
        if ctr > maxit:
            jerr = -l
    # ------------------final update & save results-------------#
    if ni > pmax:
        jerr = -10000 - l
        exit()
    if ni > 0:
        beta[1:ni] = b[m[1:ni]]
        nbeta[l] = ni
        b0[l] = al
        nalam = l
        if flmin > 1.0:
            me = find_count(abs(beta[1:ni]) > 0)
        if me > dfmax:
            exit()
def fortran_prox_sgl(gstart, gend, nvars, nobs, x, r, b, al, gamma, pf, peps, gw, step):
    s = step
    bold[gstart:gend] = b[gstart:gend]
    for g in range(gstart, gend):
        u = b[g] + s * np.multiply(x[:, g], r)/nobs
        v = abs(u) - s*al*gamma*pf[g]
        if v>0.0:
            tmp = math.copysign(v, u)
        else:
            tmp = 0.0
        b[g] = tmp
    normg = math.sqrt(np.multiply(b[gstart:gend], b[gstart:gend]))
    maxg = 0.0
    vg = s*gw*al*(1.0-gamma)/normg
    if normg is 0.0:
        vg = big
    for g in range(gstart, gend):
        scl = 1.0 -pf(g)*vg
        scl = max(scl, 0.0)
        tmp = scl*b[g]
        d = tmp - bold[g]
        r = r - x[:, g]*d
        maxg = max(maxg, abs(d))
        b[g] = tmp
    if maxg < peps:
        exit()
def computegapobs(whichfoldnot, N, l):
    out = []
    band = list(range((whichfoldnot - l), (whichfoldnot + l)))
    vec = list(range(1, N + 1))
    for i in range(len(vec)):
        index = 0
        for j in range(len(band)):
            if vec[i] == band[j]:
                index = -1
        if index != -1:
            out.append(vec[i])
    out = list(set(out))
    return out


def find_max(x):
    max_val = 0
    for i in range(len(x)):
        if x[i] > max_val:
            max_val = x[i]
    return max_val


def find_diff(x):
    diff_arr = []
    for i in range(len(x) - 1):
        diff_arr.append(x[i + 1] - x[i])
    return diff_arr


def find_index(x):
    index_arr = []
    for i in range(len(x)):
        if x[i] == 1:
            index_arr.append(i)
    index_arr.append(len(x))
    return index_arr


def find_whichgaptrain(x, y):
    del x[y[0] - 1:y[len(y) - 1]]
    return x


if __name__ == "__main__":
    random.seed(123)
    x = np.random.rand(100, 20)
    x = np.asmatrix(x)
    beta = [[5], [4], [3], [2], [1], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0], [0]]
    beta = np.asmatrix(beta)
    y = x * beta
    gindex = [1, 1, 1, 1, 1, 2, 2, 2, 2, 2, 3, 3, 3, 3, 3, 4, 4, 4, 4, 4]
    gindex = np.asmatrix(gindex)
    gindex = np.squeeze(np.asarray(gindex))
    gamma = 0.5
    standardize = False
    intercept = True
    seed = 1234
    K = 99
    nf = None
    nlambda = 100
    lamb = None

    tscv_sglfit(x, y, nlambda, 'single', nf, lamb, gindex, standardize, intercept)
