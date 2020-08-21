# -*- coding: utf-8 -*-
"""
Created on Sun Aug 16 10:20:28 2020

@author: Laptop
"""

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
  
  fit = np.array(range(63)).reshape(21,3)
  alphahat = fit[0:zparam,]
  betahat  = fit[zparam:,]
  return {'alphahat':alphahat,'betahat':betahat}

    
x = np.array(np.random.normal(0, 1, 100*20)).reshape(100,20)
y = np.random.normal(0, 1, 100)
index =list(range(1,21))
Z =  np.ones([len(y),1])
print(sgl_fit(X = x, Z = Z, y = y, index = index, lambdas =list(range(1,4)), gamma_w = None)['alphahat'] )
print(sgl_fit(X = x, Z = Z, y = y, index = index, lambdas = list(range(1,4)), gamma_w = 1)['betahat'] )

full_est=sgl_fit(X = x, Z = Z, y = y, index = index, lambdas =list(range(1,4)), gamma_w = 1)

alpha_path = full_est['alphahat']
beta_path = full_est['betahat']
i=0
tmp = lambda: None
tmp.alpha = alpha_path.reshape(1,alpha_path.shape[1]*alpha_path.shape[0],order='F')[0,i]
tmp.beta = beta_path[:,i] 

object=tmp

newX=x
 