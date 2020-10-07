#ifndef ALLFUNCTIONS_H
#define ALLFUNCTIONS_H
#include <armadillo>
using namespace arma;
#ifndef LD_ESTIMATE_SKELETON
#define LD_ESTIMATE_EKELETON
arma::vec fastols(const arma::vec & Y, const arma::mat & X, double intercept);
arma::vec fastals(const arma::vec & Y, const arma::mat & X, double intercept, double tau, double maxIter, double thresh);
arma::mat boundc(colvec x, colvec dx);
arma::vec fastrq(const arma::vec & Y, const arma::mat & X, double intercept, double tau);
arma::vec midas_pr(const arma::vec & Y, const arma::mat & X, double intercept, double tau, double which_loss, double num_evals);
arma::vec midasar_pr(const arma::vec & Y, const arma::mat & YLAG, const arma::mat & X, double intercept, double tau, double which_loss, double num_evals);
#endif
#ifndef FIT_SGL_SKELETON
#define FIT_SGL_SKELETON
void linGradCalc(int *nrow, double *eta, double *y, double *ldot);
double linNegLogLikelihoodCalc(int *nrow, double *eta, double *y);
void linSolverLamGrid(double *X, double *y, int* index, int *nrow, int *ncol, int *numGroup, double *beta, int *rangeGroupInd, int *groupLen, double *lambda1, double *lambda2, int *innerIter, double *thresh, double *ldot, double *nullBeta, double *gamma, double *eta, int* betaIsZero, int& groupChange, int* isActive, int* useGroup, double *step, int *reset);
void linNestLamGrid(double *X, double*y, int* index, int *nrow, int *ncol, 
                    int *numGroup, int *rangeGroupInd, int *groupLen, 
                    double *lambda1, double *lambda2, double *beta, 
                    int *innerIter, int *outerIter, double *thresh, 
                    double *outerThresh, double *eta, double *gamma, 
                    int *betaIsZero, double *step, int *reset);
arma::vec cvfolds(double nfolds, double nrow);
double getmin_cpp(arma::vec lambda,arma::vec cvm,arma::vec cvsd, int which_lambda);
arma::mat cpp_sgl_fit(arma::vec& beta0, arma::mat& Z, arma::mat& X, arma::vec& y, arma::vec& index,
                      arma::vec& lambda1, arma::vec& lambda2, 
                      int innerIter, int outerIter, double thresh, 
                      double outerThresh, double gamma_solver, 
                      double step, int reset);
arma::mat cpp_sgl_fitpath(arma::mat& X, arma::mat& Z, arma::vec& y, arma::vec& index, double dummies,
                          arma::vec l1_frac, arma::vec l21_frac, arma::vec dummies_index,
                          arma::vec& lambdas, double gamma_w, 
                          int innerIter, int outerIter, double thresh, 
                          double outerThresh,double gamma_solver, 
                          double step, int reset);
#endif
#endif