# -*- coding: utf-8 -*-
"""
Created on Sat Aug 15 19:07:33 2020

@author: Laptop
"""

import numpy as np
import numpy as np
import pandas as pd
from sty import fg, bg, ef, rs
import math
import warnings
import sys
from statistics import mean






def sgl_fit(X, Z, y, index, lambdas, gamma_w=None, l1_factor=None, l21_factor=None,
            dummies_index=None, inner_iter=None, outer_iter=None,
            thresh=None, outer_thresh=None):
  
  X = np.array(X)
  Z = np.array(Z)
  y = y #check
  num_groups = max(index)
  if(lambdas==None):
    sys.exit('[Error] lambda sequence was not specified')
  
  if(pd.isnull(l1_factor)):
    l1_factor = np.zeros([1,Z.shape[1]])
  if(pd.isnull(l21_factor)):
    l21_factor = np.zeros([1,Z.shape[1]])  
  if(pd.isnull(dummies_index)): 
    dummies_index = 1
  
  if(pd.isnull(gamma_w)):
    print(fg.red + 'sg-lasso relative weight was not set. setting to default value 0.8'  +  fg.rs)
    gamma_w = 0.8
  
  
  if(pd.isnull(inner_iter)):
    inner_iter = 1e2
  
  if(pd.isnull(outer_iter)):
    outer_iter = 1e4
  
  if(pd.isnull(thresh)):
    thresh = 1e-2
  
  if(pd.isnull(outer_thresh)):
    outer_thresh = 1e-4
  
  gamma_solver = 0.8
  step = 1
  reset = 10
  
  zparam = Z.shape[1]
  xparam = X.shape[1]
  dummies = 1
  #check
  #------------------ sgl fit -----------#
  '''
  fit = cpp_sgl_fitpath(np.array(X),np.array(Z),y, index, float(dummies), 
                 l1_factor, l21_factor,dummies_index,
                  lambdas, float(gamma_w), inner_iter, int(outer_iter), float(thresh), float(outer_thresh), 
                  float(gamma_solver), float(step), int(reset))
  '''
  
  ''' THE cpp_sgl_fitpath FUNCTION IS IMPLEMENTED IN C, THEREFORE I CANNOT PROGRAM IT'''
  
  fit = np.array(range(63)).reshape(21,3)
  alphahat = fit[0:zparam,]
  betahat  = fit[zparam:,]
  return {'alphahat':alphahat,'betahat':betahat}











def path_calc(X, Z, y, index, gamma_w=None, l1_factor=0, l21_factor=0, dummies_index=1, nlam=None, min_frac=None):
  if(gamma_w==None):
    print(fg.red + 'sg-lasso relative weight was not set. setting to default value 0.8' + fg.rs)
    gamma_w = 0.8
  
  if(nlam==None):
    nlam = 100
  
  if (min_frac==None):
    k = X.shape[1]
    t = X.shape[0]
    if (k>t):
      min_frac = 0.0001
    else:
      min_frac = 0.001
    
  
  X = np.array(X)
  Z = np.array(Z)
 # y = y
  n = X.shape[0]
  max_lam_tmp = max(np.dot( np.concatenate((Z,X),axis=1).T, y ))/n
 
  # we just need approximate solutions so coordinate-descent parameters are set to be loose
  inner_iter = 1e3
  outer_iter = 1e3
  thresh = 1e-3
  outer_thresh = 1e-3
  
  # fit the initial solution
  path_est = sgl_fit(X, Z, y, index, max_lam_tmp, gamma_w, l1_factor, l21_factor,
                     dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, 
                     thresh=thresh, outer_thresh=outer_thresh)['betahat']
  tmp_lambda = max_lam_tmp
  notfound = True
  k = 0
  if np.any(path_est!=0):
    factor = 1.01
    while (notfound):
      k = k + 1
      tmp_lambda = tmp_lambda*factor
      tmp_fit = sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor,
                        dummies_index, inner_iter=inner_iter, outer_iter=outer_iter, 
                        thresh=thresh, outer_thresh=outer_thresh)['betahat']
      notfound = np.any(tmp_fit!=0)
      if (k > 1e3):
        notfound = False
        tmp_lambda = max_lam_tmp
        
        print(fg.red + 'Warning message: \n In path_calc(X, Z, y, index, gamma_w = NULL, l1_factor = 0, l21_factor = 0,  :\n  more than 1e3 iterations were used to fine-tune lambdas sequence. initial value set to max(Xty)/n' + fg.rs)
        
  else:
    factor = 0.99
    while (notfound):
      k = k + 1
      tmp_lambda = tmp_lambda*factor
      tmp_fit = sgl_fit(X, Z, y, index, tmp_lambda, gamma_w, l1_factor, l21_factor,
                        dummies_index, inner_iter=inner_iter, outer_iter=outer_iter,
                        thresh=thresh, outer_thresh=outer_thresh)['betahat']
      notfound = np.all(tmp_fit==0)
      if (k > 1e3):
        notfound = False
        tmp_lambda = max_lam_tmp
        print(fg.red + 'Warning message: \n In path_calc(X, Z, y, index, gamma_w = NULL, l1_factor = 0, l21_factor = 0,  :\n  more than 1e3 iterations were used to fine-tune lambdas sequence. initial value set to max(Xty)/n' + fg.rs)
    tmp_lambda = tmp_lambda/0.99
  max_lam = tmp_lambda
  min_lam = min_frac*max_lam
  log_min_lam = math.log(min_lam)
  log_max_lam = math.log(max_lam)
  seq_loglambdas = list(np.linspace(log_max_lam,log_min_lam, num = nlam)) 
  lambdas = np.exp(seq_loglambdas)
  return lambdas







def predict_reg_sgl(object, newX):
  alpha = object['alpha']
  beta = object['beta']
  predZ = alpha
  predX = list(np.dot(newX,beta))
  pred = list(predZ+predX)
  return {'pred': pred, 'predZ': predZ, 'predX' : predX}

def reg_sgl(X, y, index, gamma_w=None, full_est=None, method_choice=("ic","cv","initial"), 
             nlam=100, lambdas=None,min_frac=None, nfolds=10, lambda_choice=("min","1se"),
             ic_choice=("bic","aic","aicc"),
             num_cores = None, verbose=False,thresh=None, outer_thresh=None, inner_iter=None,
             outer_iter=None):
      # check if input data has no na entries
    if(pd.isnull(y).any()):
        sys.exit('[Error] y has NA entries, check and rerun')
   
    if(pd.isnull(X).any()):
        sys.exit('[Error] X has NA entries, check and rerun')
   
     # get settings
    method = method_choice
    if(method not in ("ic","cv","initial")):
      sys.exit('[Error] Choose only one of these 3: "ic","cv","initial"')
      
    lambda_choice = lambda_choice
    if(lambda_choice not in ("min","1se")):
      sys.exit('[Error] Choose only one of these 2: "min","1se"')    
    
    ic_choice = ic_choice
    if(ic_choice not in ("bic","aic","aicc")):
      sys.exit('[Error] Choose only one of these 3: "bic","aic","aicc"')    
        
    
    if(lambda_choice=="min"):
        which_lambda = float(0)
    if (lambda_choice=="1se"):
        which_lambda = float(1)
     
    tp =X.shape
    t =tp[0]
    p =tp[1]
     
    # intercept dummy
    Z =np.ones([len(y),1]) 
  # zero penalty weight for the intercept
    l1_factor = 0
    l21_factor = 0 
    dummies_index = 1 # group index for the intercept
    X =  np.array(X)
   # y =y #check
    if (pd.isnull(lambdas)):
        if(verbose):
             print(fg.red + 'computing lambda sequence' + fg.rs)
        lambdas =path_calc(X, Z, y, index, gamma_w, l1_factor, l21_factor,
                           dummies_index)
  # compute entire solution
    if (pd.isnull(full_est)):
        if(verbose):
            print(fg.red + 'computing full path solution' + fg.rs)
        full_est =sgl_fit(X, Z, y, index, lambdas, gamma_w, l1_factor, l21_factor,
                         dummies_index, inner_iter, outer_iter, thresh, outer_thresh)
  
    if (method=="initial"):
        if(verbose):
            print(fg.red + 'returning initial full path estimates and lambda sequence' + fg.rs)
        return {'full_est':full_est, 'lambdas':lambdas}
    
    if (method=="ic"):
        if(verbose):
            print(fg.red + 'computing solution using - ' + ic_choice + ' - as information criterion' +  fg.rs)
        crit = np.zeros(nlam)
        alpha_path = full_est['alphahat']
        beta_path = full_est['betahat']
  
        for i in range(0,nlam):
            tmp = {}
            tmp['alpha'] = alpha_path.reshape(1,alpha_path.shape[1]*alpha_path.shape[0],order='F')[0,i]
            tmp['beta'] = beta_path[:,i] 
      # compute df as |beta|_0
            beta_0_card = np.sum(tmp.beta != 0)
            df = beta_0_card
      # compute penalty
            if (ic_choice=="bic"): pen =  math.log(t)/t*df #check
            if (ic_choice=="aic") : pen = 2/t*df
            if (ic_choice=="aicc"):  pen = 2*df/(t - df - 1)
            yhat = predict_reg_sgl(object = tmp, newX = X)['pred']
            sigsqhat =sum(np.power(y-mean(y),2))/t
            mse = sum(np.power(y-yhat,2))/t    
            crit[i] =mse/sigsqhat + pen 
            
        lambda_ic = lambdas[np.where(crit==min(crit))]
        
    #take single solution in case many are optimum:
        if (len(lambda_ic)>1):
            lambda_ic =min(lambda_ic)
        min_crit =list(crit)
        alpha =full_est['alphahat'].reshape(1,full_est['alphahat'].shape[1]*full_est['alphahat'].shape[0],order='F')[0,np.where(lambdas==lambda_ic)]
        beta =full_est['betahat'][:,np.where(lambdas==lambda_ic)]
       
    return{'alpha_path' : alpha_path, 'beta_path' : beta_path, 'lambdas' : lambdas, 
          'alpha' : alpha, 'beta' : beta, 'lambda_ic' : lambda_ic,
          'min_crit' : min_crit, 'full_est' : full_est}
   
    if (method=="cv") :
        if(verbose):
            print(fg.red + 'computing' + nfolds + '-fold cross-validation' +  fg.rs)
    # compute folds
    '''
       is in C
       foldsid =cvfolds(float(nfolds), float(t))+1'''
       
    foldsid = np.array([1,1,1,1,1,1,1,1,1,1,2,2,2,2,2,2,2,2,2,2,3,3,3,3,3,3,3,3,3,3,4,4,4,4,4,4,4,4,4,4,5,5,5,5,5,5,5,5,5,5,
                6,6,6,6,6,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,10,10]).reshape(100,1)
     #======== mc stuff ========# 
    if (pd.isnull(num_cores)):
        from platform import system
        import multiprocessing
        
        if(system()=="Unix" or system()=="unix"):
          num_cores =multiprocessing.cpu_count()
        if(system()=="Windows"):
          num_cores =multiprocessing.cpu_count()
    
    #cl =parallel::makeCluster(num_cores)
    #doSNOW::registerDoSNOW(cl)
     if(verbose) :
     # pb =utils::txtProgressBar(max=nfolds, style=3)
     # progress =function(n) utils::setTxtProgressBar(pb, n)
     # opts =list(progress=progress)
       else :
      opts =None
    
    
    #======== Main loop ========# 
   # output =foreach::foreach(k = 1:nfolds, .packages = c("midasml"), .options.snow=opts) %dopar% {
        # computing CV solutions
       # which_val =np.where(foldsid==k)
        #which_train <-np.where(foldsid!=k)

        #y_train =y[which_train]
        #X_train =X[which_train, ]
        #Z_train =Z[which_train, ]
        
        #est_new =sgl_fit(X_train, Z_train, y_train, index, lambdas, gamma_w, l1_factor, l21_factor, dummies_index)
        #est_new
    
    parallel::stopCluster(cl)
    
    # initialize cv mean and sd storage
    cvm =pd.zeros(nlam)
    #cvsd_fold =matrix(0, nrow=nlam, ncol=nfolds)
    for (k in 1:nfolds):
      #tmp = output[[k]]
      est_new = np.concatenate((tmp['alphahat'],tmp['betahat']), axis=0)    #rbind(tmp$alphahat,tmp$betahat)
      which_val = np.where(foldsid==k)
      which_train = np.where(foldsid!=k)
      
      y_val = y[which_val]
      X_val = X[which_val, ]
      Z_val = Z[which_val, ]
      ZX_val = as.matrix(cbind(Z_val,X_val))
      tmp = updatecvmsd(as.vector(cvm),  as.matrix(cvsd_fold), as.double(nlam), as.matrix(est_new), as.double(k), as.vector(y_val), as.matrix(ZX_val))
      cvm = tmp$cvm
      cvsd = tmp$cvsd_fold
    
    
    cvm = cvm/nfolds
    cvsd = apply(cvsd,1,sd) * sqrt(nfolds)
    lambda_cv = getmin_cpp(lambdas, cvm, cvsd, which_lambda)
    min_crit = list(cvm = cvm, cvsd = cvsd)
    alpha_path = full_est$alphahat
    beta_path = full_est$betahat
    alpha = full_est$alphahat[np.where(lambdas==lambda_cv)]
    beta = full_est$betahat[,np.where(lambdas==lambda_cv)]
    return {'alpha_path': alpha_path, 'beta_path': beta_path', 'lambdas': lambdas, 
            'alpha': alpha, 'beta': beta, 'lambda_cv': lambda_cv, 
            'min_crit': min_crit, 'full_est': full_est}
    