include "cyarma.pyx"
cimport numpy as np
import numpy as np
#Code is going here
#External functions from allfunctions.h
cdef extern from "allfunctions.h":
  void linGradCalc(int* nrow, double *eta, double *y, double *ldot)
  double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y)
  void linSolverLamGrid(double *X, double *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset)
  void linNestLamGrid(double *X, double*y, int* index, int *nrow, int *ncol,
                    int *numGroup, int *rangeGroupInd, int *groupLen,
                    double *lambda1, double *lambda2, double *beta,
                    int *innerIter, int *outerIter, double *thresh,
                    double *outerThresh, double *eta, double *gamma,
                    int *betaIsZero, double *step, int *reset)
  vec cvfolds(double nfolds, double nrow)
  double getmin_cpp(vec lambda_,vec cvm,vec cvsd, int which_lambda);
  mat cpp_sgl_fit(vec& beta0, mat& Z, mat& X, vec& y, vec& index,
                      vec& lambda1, vec& lambda2,
                      int innerIter, int outerIter, double thresh,
                      double outerThresh, double gamma_solver,
                      double step, int reset)
  mat cpp_sgl_fitpath(mat& X, mat& Z, vec& y, vec& index, double dummies,
                          vec l1_frac, vec l21_frac, vec dummies_index,
                          vec& lambdas, double gamma_w,
                          int innerIter, int outerIter, double thresh,
                          double outerThresh,double gamma_solver,
                          double step, int reset)
  vec fastols(const vec & Y, const mat & X, double intercept)
  vec fastals(const vec & Y, const mat & X, double intercept, double tau, double maxIter, double thresh)
 # mat boundc(colvec x, colvec dx) #not wrapped in python some api issue
  vec fastrq(const vec & Y, const mat & X, double intercept, double tau)
  vec midas_pr(const vec & Y, const mat & X, double intercept, double tau, double which_loss, double num_evals)
  vec midasar_pr(const vec & Y, const mat & YLAG, const mat & X, double intercept, double tau, double which_loss, double num_evals)
#Building wrapper for external functions with support of numpy
def lin_grad_calc(np.ndarray[np.int_t,ndim=1] nrow,np.ndarray[np.double_t,ndim=1] eta,np.ndarray[np.double_t,ndim=1] y,np.ndarray[np.double_t,ndim=1] ldot):
  linGradCalc(<int*>nrow.data,<double*>eta.data,<double*>y.data,<double*>ldot.data)

def lin_neg_log_likelihood_calc(np.ndarray[np.int_t,ndim=1] nrow,np.ndarray[np.double_t,ndim=1] eta,np.ndarray[np.double_t,ndim=1] y):
  return linNegLogLikelihoodCalc(<int*>nrow.data,<double*>eta.data,<double*>y.data)

def lin_solver_lam_grid(
        np.ndarray[np.double_t,ndim=1] X,np.ndarray[np.double_t,ndim=1] y,
        np.ndarray[np.int_t,ndim=1] index,np.ndarray[np.int_t,ndim=1] nrow,
        np.ndarray[np.int_t,ndim=1] ncol,np.ndarray[np.int_t,ndim=1] numGroup,
        np.ndarray[np.double_t,ndim=1] beta,np.ndarray[np.int_t,ndim=1] rangeGroupInd,
        np.ndarray[np.int_t,ndim=1] groupLen,np.ndarray[np.double_t,ndim=1] lambda1,np.ndarray[np.double_t,ndim=1] lambda2,
        np.ndarray[np.int_t,ndim=1] innerIter,np.ndarray[np.double_t,ndim=1] thresh,np.ndarray[np.double_t,ndim=1] ldot,
        np.ndarray[np.double_t,ndim=1] nullBeta,np.ndarray[np.double_t,ndim=1] gamma,
        np.ndarray[np.double_t,ndim=1] eta,np.ndarray[np.int_t,ndim=1] betaIsZero,groupChange,np.ndarray[np.int_t,ndim=1] isActive,
        np.ndarray[np.int_t,ndim=1] useGroup,np.ndarray[np.double_t,ndim=1] step,np.ndarray[np.int_t,ndim=1] reset
):
    linSolverLamGrid(
        <double*>X.data,<double*>y.data,<int*>index.data,
        <int*>nrow.data,<int*>ncol.data,<int*>numGroup.data,
        <double*>beta.data,<int*>rangeGroupInd.data,<int*>groupLen.data,
        <double*>lambda1.data,<double*>lambda2.data,<int*>innerIter.data,
        <double*>thresh.data,<double*>ldot.data,<double*>nullBeta.data,
        <double*>gamma.data,<double*>eta.data,<int*>betaIsZero.data,
        groupChange,<int*>isActive.data,<int*>useGroup.data,<double*>step.data,
        <int*>reset.data
    )

def lin_nest_lam_grid(
   np.ndarray[np.double_t,ndim=1] X,np.ndarray[np.double_t,ndim=1] y,
   np.ndarray[np.int_t,ndim=1] index,np.ndarray[np.int_t,ndim=1] nrow,
   np.ndarray[np.int_t,ndim=1] ncol,np.ndarray[np.int_t,ndim=1] numGroup,
   np.ndarray[np.int_t,ndim=1] rangeGroupInd,np.ndarray[np.int_t,ndim=1] groupLen,
   np.ndarray[np.double_t,ndim=1] lambda1,np.ndarray[np.double_t,ndim=1] lambda2,
   np.ndarray[np.double_t,ndim=1] beta,np.ndarray[np.int_t,ndim=1] innerIter,
   np.ndarray[np.int_t,ndim=1] outerIter,np.ndarray[np.double_t,ndim=1] thresh,
   np.ndarray[np.double_t,ndim=1] outerThresh,np.ndarray[np.double_t,ndim=1] eta,
   np.ndarray[np.double_t,ndim=1] gamma,np.ndarray[np.int_t,ndim=1] betaIsZero,
   np.ndarray[np.double_t,ndim=1] step,np.ndarray[np.int_t,ndim=1] reset
):
  linNestLamGrid(
      <double*>X.data,<double*>y.data,
      <int*>index.data,<int*>nrow.data,
      <int*>ncol.data,<int*>numGroup.data,
      <int*>rangeGroupInd.data,<int*>groupLen.data,
      <double*>lambda1.data,<double*>lambda2.data,
      <double*>beta.data,<int*>innerIter.data,
      <int*>outerIter.data,<double*>thresh.data,
      <double*>outerThresh.data,<double*>eta.data,
      <double*>gamma.data,<int*>betaIsZero.data,
      <double*>step.data,<int*>reset.data
  )

def cvfold_s(nfolds,nrow):
    cdef vec result=cvfolds(nfolds,nrow)
    cdef double* Xptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,dtype=np.double)
    cdef double* Dptr=<double*>D.data
    for i in range(result.n_elem):
        Dptr[i]=Xptr[i]
    return D

def getmin(np.ndarray[np.double_t,ndim=1] lambda_,np.ndarray[np.double_t,ndim=1] cvm,np.ndarray[np.double_t,ndim=1] cvsd,which_lambda):
    cdef vec* lambda_Vec=new vec(<double*>lambda_.data,lambda_.shape[0],False,True)
    cdef vec* cvm_Vec=new vec(<double*>cvm.data,cvm.shape[0],False,True)
    cdef vec* cvsd_Vec=new vec(<double*>cvsd.data,cvsd.shape[0],False,True)
    return getmin_cpp(<vec>lambda_Vec[0],<vec>cvm_Vec[0],<vec>cvsd_Vec[0],which_lambda)

def sgl_fit(
        np.ndarray[np.double_t,ndim=1] beta0,np.ndarray[np.double_t,ndim=2] Z,
        np.ndarray[np.double_t,ndim=2] X,np.ndarray[np.double_t,ndim=1] y,
        np.ndarray[np.double_t,ndim=1] index,np.ndarray[np.double_t,ndim=1] lambda1,
        np.ndarray[np.double_t,ndim=1] lambda2,innerIter,outerIter,thresh,
        outerThresh,gamma_solver,step,reset
):
  cdef vec* beta0_Vec=new vec(<double*>beta0.data,beta0.shape[0],False,True)
  if not Z.flags.f_contiguous:
      Z=Z.copy(order='F')
  if not X.flags.f_contiguous:
      X=X.copy(order='F')
  cdef mat* ZMat=new mat(<double*>Z.data,Z.shape[0],Z.shape[1],False,True)
  cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)
  cdef vec* y_Vec=new vec(<double*>y.data,y.shape[0],False,True)
  cdef vec* index_Vec=new vec(<double*>index.data,index.shape[0],False,True)
  cdef vec* lambda1_Vec=new vec(<double*>lambda1.data,lambda1.shape[0],False,True)
  cdef vec* lambda2_Vec=new vec(<double*>lambda2.data,lambda2.shape[0],False,True)

  cdef mat result=cpp_sgl_fit(
    beta0_Vec[0],
    ZMat[0],
    XMat[0],
    y_Vec[0],
    index_Vec[0],
    lambda1_Vec[0],
    lambda2_Vec[0],
    innerIter,outerIter,thresh,
    outerThresh,gamma_solver,step,reset)
  cdef np.ndarray[np.double_t,ndim=2] Dx=np.empty((result.n_rows,result.n_cols),dtype=np.double)
  cdef double* mptr=result.memptr()
  cdef double* pxptr=<double*> Dx.data
  for i in range(result.n_rows*result.n_cols):
      pxptr[i]=mptr[i]
  return Dx

def sgl_fitpath(
    np.ndarray[np.double_t,ndim=2] X,np.ndarray[np.double_t,ndim=2] Z,
    np.ndarray[np.double_t,ndim=1] y,np.ndarray[np.double_t,ndim=1] index,dummies,
    np.ndarray[np.double_t,ndim=1] l1_frac,np.ndarray[np.double_t,ndim=1] l21_frac,
    np.ndarray[np.double_t,ndim=1] dummies_index,np.ndarray[np.double_t,ndim=1] lambdas,
    gamma_w,innerIter,outerIter,thresh,outerThresh,gamma_solver,step,reset
):
    if not X.flags.f_contiguous:
        X=X.copy(order='f')
    if not Z.flags.f_contiguous:
        Z=Z.copy(order='f')
    
    cdef mat*  XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)
    cdef mat*  ZMat=new mat(<double*>Z.data,Z.shape[0],Z.shape[1],False,True)
    cdef vec* YVec=new vec(<double*>y.data,y.shape[0],False,True)
    cdef vec* indexVec=new vec(<double*>index.data,index.shape[0],False,True)
    cdef vec* l1_frac_Vec=new vec(<double*>l1_frac.data,l1_frac.shape[0],False,True)
    cdef vec* l21_frac_Vec=new vec(<double*>l21_frac.data,l21_frac.shape[0],False,True)
    cdef vec* lambdas_Vec=new vec(<double*>lambdas.data,lambdas.shape[0],False,True)
    cdef vec* dummies_indexVec=new vec(<double*>dummies_index.data,dummies_index.shape[0],False,True)
    
    cdef mat result= cpp_sgl_fitpath(
          XMat[0],
          ZMat[0],
          YVec[0],
          indexVec[0],
          dummies,
          l1_frac_Vec[0],
          l21_frac_Vec[0],
          dummies_indexVec[0],
          lambdas_Vec[0],
          gamma_w,innerIter,outerIter,thresh,outerThresh,gamma_solver,step,reset
      )
    cdef double* rptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=2] D=np.empty((result.n_rows,result.n_cols),dtype=np.double)
    cdef double* DPtr=<double*>D.data
    for i in range(result.n_cols*result.n_rows):
        DPtr[i]=rptr[i]
    return D

def fostol_s(
        np.ndarray[np.double_t,ndim=1] Y,
        np.ndarray[np.double_t,ndim=2] X,
        intercept
):
    cdef vec* Y_Vec=new vec(<double*>Y.data,Y.shape[0],False,True)
    if not X.flags.f_contiguous:
        X=X.copy(order='F')
    cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)

    cdef vec result=fastols(<vec>Y_Vec[0],<mat>XMat[0],intercept)
    cdef double* mptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,dtype=np.double)
    cdef double* DXPTR=<double*>D.data
    for i in range(result.n_elem):
        DXPTR[i]=mptr[i]
    return D

def fastal_s(
        np.ndarray[np.double_t,ndim=1] Y,
        np.ndarray[np.double_t,ndim=2] X,
        intercept,tau,maxIter,thresh
):
    cdef vec* Y_Vec=new vec(<double*>Y.data,Y.shape[0],False,True)
    if not X.flags.f_contiguous:
        X=X.copy(order='F')
    cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)

    cdef vec result=fastals(
        Y_Vec[0],
        XMat[0],
        intercept,tau,maxIter,thresh
    )
    cdef double* mptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,np.double)
    cdef double* DXPTR=<double*>D.data
    for i in range(result.n_elem):
        DXPTR[i]=mptr[i]
    return D

def fastr_q(
        np.ndarray[np.double_t,ndim=1] Y,
        np.ndarray[np.double_t,ndim=2] X,
        intercept,tau
):
    cdef vec* Y_Vec=new vec(<double*>Y.data,Y.shape[0],False,True)
    if not  X.flags.f_contiguous:
        X=X.copy(order='F')
    cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)

    cdef vec result=fastrq(
        Y_Vec[0],
        XMat[0],
        intercept,tau
    )
    cdef double* mptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,dtype=np.double)
    cdef double* DXPTR=<double*>D.data
    for i in range(result.n_elem):
        DXPTR[i]=mptr[i]
    return D

def midaspr(
        np.ndarray[np.double_t,ndim=1] Y,
        np.ndarray[np.double_t,ndim=2] X,
        intercept,tau,which_loss,num_evals
):
    cdef vec* Y_Vec=new vec(<double*>Y.data,Y.shape[0],False,True)
    if not X.flags.f_contiguous:
        X=X.copy(order='F')
    cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)

    cdef vec result=midas_pr(
        Y_Vec[0],
        XMat[0],
        intercept,tau,which_loss,num_evals
    )
    cdef double* mptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,dtype=np.double)
    cdef double* xptr=<double*>D.data
    for i in range(result.n_elem):
        xptr[i]=mptr[i]
    return D
def midasarpr(
        np.ndarray[np.double_t,ndim=1] Y,
        np.ndarray[np.double_t,ndim=2] YLAG,
        np.ndarray[np.double_t,ndim=2] X,
        intercept,tau,which_loss,num_evals
):
    cdef vec* Y_Vec=new vec(<double*>Y.data,Y.shape[0],False,True)
    if not YLAG.flags.f_contiguous:
        YLAG=YLAG.copy(order='F')
    if not X.flags.f_contiguous:
        X=X.copy(order='F')
    cdef mat* YLAGMAT=new mat(<double*>YLAG.data,YLAG.shape[0],YLAG.shape[1],False,True)
    cdef mat* XMat=new mat(<double*>X.data,X.shape[0],X.shape[1],False,True)
    cdef vec result=midasar_pr(
            Y_Vec[0],
            YLAGMAT[0],
            XMat[0],
            intercept,tau,which_loss,num_evals
        ),
    cdef double* mptr=result.memptr()
    cdef np.ndarray[np.double_t,ndim=1] D=np.empty(result.n_elem,dtype=np.double)
    cdef double* DXPTR=<double*>D.data
    for i in range(result.n_elem):
        DXPTR[i]=mptr[i]
    return D
#Code is ended here
