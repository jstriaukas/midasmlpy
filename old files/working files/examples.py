import midasml.midasml as ml
import numpy as np

'''
ml.lin_grad_calc example take _nrow ,_eta, _y,and _ldot
calculate and write result into _ldot back same as void linGradCalc in c++
'''
_ldot=np.array([0,0,0,0,0],dtype='float64')
ml.lin_grad_calc(np.array([5]),np.array([2,5,4,6,8],dtype='float64'),np.array([6,5,8,9,5],dtype='float64'),_ldot)
print(_ldot)
'''
ml.lin_neg_log_likelihood_calc example take _nrow,_eta, and _y and return double type calculated value
same as double linNegLogLikelihoodCalc(...) in c++
'''
value=ml.lin_neg_log_likelihood_calc(np.array([5]),_ldot,np.array([1,2,3,6,8],dtype='float64'))
print(value)
'''
Do it your self because i don't know how much value provide for computation
'''

_X=np.array([1,2,5,4,4,8],dtype='float64')#provide list of value in [],
_y=np.array([1,2,5,4,5,6],dtype='float64')#provide here
_index=np.array([0,1,2,3,4,5],dtype='int')# provide here
_nrow=np.array([3],dtype='int')# provide one value inside [] e.g. [56]
_ncol=np.array([3],dtype='int')#  provide one value inside []
_numGroup=np.array([1],dtype='int') # provide value here
_beta=np.array([0,0,0],dtype='float64') # provide here list of value inside []
_rangeGroupInd=np.array([2],dtype='int') # provide value here
_groupLen=np.array([1],dtype='int')# provide value
_lambda1=np.array([0,0,0],dtype='float64')# provide value
_lambda2=np.array([1,1,1],dtype='float64')# provide value
_innerIter=np.array([10],dtype='int')#provide value
_thresh=np.array([5],dtype='float64') # provide value here
_ldot=np.array([1],dtype='float64') # provide value
_nullBeta=np.array([0,1,2],dtype='float64') # provide value
_gamma=np.array([2,2,5,5],dtype='float64') # provide value
_eta=np.array([5,4,8,7],dtype='float64') # provide value
_betaIsZero=np.array([5,5,4,5],dtype='int') # provide value
_groupChange=5 # provide a double value or may be it can be a output
_isActive=np.array([1],dtype='int')# provide value
_useGroup=np.array([1],dtype='int')# provide value
_step=np.array([3],dtype='float64')# provide value
_reset=np.array([1],dtype='int')#provide value
ml.lin_solver_lam_grid(
    _X,_y,_index,_nrow,_ncol,_numGroup,_beta,_rangeGroupInd,
    _groupLen,_lambda1,_lambda2,_innerIter,_thresh,_ldot,_nullBeta,
    _gamma,_eta,_betaIsZero,_groupChange,_isActive,_useGroup,_step,_reset
)
print(_ldot)
print(_groupChange)
print(ml.cvfolds(500,1000))
