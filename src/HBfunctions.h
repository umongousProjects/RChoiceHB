#ifndef HBFUNCTIONS_H
#define HBFUNCTIONS_H

#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h> //M_PI
#include <algorithm> // std::unique, std::distance
#include <tuple>//std::tuple
#include <chrono>//high_resolution_clock
#include <vector>


using namespace Rcpp;
using namespace arma;
using namespace std::chrono;

void coutMat(mat &thisM,char  text[]);

void rGibbs(const mat &rGibbs_y,const  int nu,const  mat V,const  int nPar,const  int nRespondents,mat &rmultireg_B, 
                mat &rooti, Environment &myEnv,double IRd);

mat doV(const int nPar,const  int numOfVars,const  IntegerVector &featureLengths,const  int nu);

mat doX(const int nPar,const  int numOfVars,const  IntegerVector &featureLengths,const  mat &X_cat);

mat getLLMnl(const mat &X,const  int nPar,const  mat &newbeta,const  int packagesPerScreen,const  uvec &yIndices,
        const  int numberOfScreens);

mat getLndMvn(const mat &x,const mat &mu,const mat &rooti,const  int nPar);

void mnlRwMetropOnce(const uvec &yIndices,const mat &X, mat &oldbeta,mat &oldll, 
        double s,const  mat &incroot, const  mat &betabar,const  mat &rootpi,const  int nRespondents,const  int numberOfScreens, 
        const int packagesPerScreen, const uvec &yIndices_old,const int nPar, int &stay);

#endif
