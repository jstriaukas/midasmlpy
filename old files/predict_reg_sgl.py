

import numpy as np
import numpy as np
import pandas as pd
from sty import fg, bg, ef, rs
import math
import warnings
import sys



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
  
  fit = np.array(np.random.normal(0, 1, 63)).reshape(21,3)
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



x = np.array(np.random.normal(0, 1, 100*20)).reshape(100,20)
y = np.random.normal(0, 1, 100)
index =list(range(1,21))
Z =  np.ones([len(y),1])

full_est=sgl_fit(X = x, Z = Z, y = y, index = index, lambdas =list(range(1,4)), gamma_w = 1)

alpha_path = full_est['alphahat']
beta_path = full_est['betahat']
i=0
tmp = {}
tmp['alpha'] = alpha_path.reshape(1,alpha_path.shape[1]*alpha_path.shape[0],order='F')[0,i]
tmp['beta'] = beta_path[:,i] 
print(predict_reg_sgl(tmp, x)['pred'])
print("\n")
print(predict_reg_sgl(tmp, x)['predZ'])
print("\n")
print(predict_reg_sgl(tmp, x)['predX'])
#yhat=predict_reg_sgl(tmp, x)['pred']



