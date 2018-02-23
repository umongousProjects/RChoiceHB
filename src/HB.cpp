#include <RcppArmadillo.h>
// [[Rcpp::depends(RcppArmadillo)]]
#include <math.h> //M_PI
#include <algorithm> // std::unique, std::distance
#include <tuple>//std::tuple
#include <chrono>//high_resolution_clock
#include <vector>

#include "HBfunctions.h"

//using namespace Rcpp;
//using namespace arma;
//using namespace std::chrono;

// [[Rcpp::export]]
arma::mat testrcpp(arma::mat XX) {
        return(XX*XX);
}

// [[Rcpp::export]]
Rcpp::List RChoiceHB(int nRespondents, int nPar,Rcpp::IntegerVector featureLengths, arma::mat X_cat, arma::uvec Y_indexed,arma::uvec yIndices_old,int numberOfScreens, int packagesPerScreen,double IRd,int rPid) {
  int numOfVars = X_cat.n_cols;
  int nu = nPar*2;
  int keep=5;
  int stay=0;
  double s = 0.1;
  
  arma::mat V = doV(nPar, numOfVars, featureLengths, nu);
  arma::mat oldBetas = arma::mat(nRespondents, nPar, fill::zeros);
  arma::mat X = doX(nPar, numOfVars, featureLengths, X_cat);
  arma::mat oldll;
  arma::mat rootPi;
  arma::mat incRoot;
  arma::mat betaBar=arma::mat(1, nPar, fill::zeros);
  arma::mat oldBetasSQ;
  double RMS;
  double RMSa;
  double acceptr;
  
  arma::mat thisbetaBar;
  arma::mat thisrootPi;
  arma::mat thisoldbetas;
  arma::mat thisoldll;
  double thisStay=0;
  
  bool converged=false;
  bool removeLastEl=false;
  std::vector<double> RMSs;
  std::vector<double> stays;
  int halfIndice;
  double RMSsMean1stHalf;
  double RMSsMean2ndHalf;
  double difference;
  int convergeStartedAt=0;
  
  arma::mat oldbetasSum;
  double oldbetasN=1.0;
  bool firstTime=true;
  
  Environment myEnv = Environment::global_env();
  Function chol2inv3 = myEnv["chol2inv3"];
  for (int iteration = 0; iteration < 100000; iteration++) {
    rGibbs(oldBetas, nu, V, nPar, nRespondents,betaBar,rootPi,myEnv,IRd);//betaBar,rootPi out
    if (iteration == 0) {
      oldll = getLLMnl(X, nPar, oldBetas, packagesPerScreen, Y_indexed, numberOfScreens);
    }
    incRoot = as < arma::mat > (chol2inv3(Rcpp::Named("x", rootPi * rootPi.t())));
    mnlRwMetropOnce(Y_indexed, X, oldBetas,oldll, s, incRoot,betaBar,rootPi, nRespondents, 
                    numberOfScreens, packagesPerScreen, yIndices_old,nPar,stay);//oldBetas,oldll,stay out
    oldBetasSQ=oldBetas;
    oldBetasSQ.for_each([](arma::mat::elem_type & val) {
      val = pow(val,2.0);
    });
    RMS=sqrt(mean(mean(oldBetasSQ)));
    if(iteration==0){
      RMSa = RMS;
    }else{
      RMSa = 0.99 * RMSa + 0.01 * RMS;
    }
    acceptr=nRespondents-stay;
    stays.push_back(stay);
    if (acceptr/nRespondents < 0.3) {
      s *= 0.9;
    } else if (acceptr/nRespondents > 0.3) {
      s *= 1.1;
    }
    if(!converged){
      if(remainder(iteration,100)==0){
        
        cout<<"session_id:"<<rPid<<" burn-in iteration:"<<iteration<<"\n"<<std::flush;
        RMSs.push_back(RMSa);
        if(removeLastEl){
          RMSs.erase(RMSs.begin());
          removeLastEl=false;
        }else{
          removeLastEl=true;
        }
        if(iteration>5000){
          halfIndice=RMSs.size()/2;
          RMSsMean1stHalf=std::accumulate(RMSs.begin(), RMSs.begin()+halfIndice, 0.0);
          RMSsMean2ndHalf=std::accumulate(RMSs.begin()+halfIndice+1, RMSs.end(), 0.0);
          RMSsMean1stHalf=RMSsMean1stHalf/halfIndice;
          RMSsMean2ndHalf=RMSsMean2ndHalf/halfIndice;
          difference=fabs(RMSsMean1stHalf-RMSsMean2ndHalf);
          if(difference<0.03 && converged==false){
            converged=true;
            convergeStartedAt=iteration;
          }
          
        }
        
      }
    }
    if (converged) {
      if(remainder(iteration,100)==0){
        cout<<"session_id:"<<rPid<<" sampling iteration:"<<iteration<<"\n"<<std::flush;
      }
      
      if(iteration-convergeStartedAt>5000){
        oldbetasSum/=oldbetasN;
        break;
        
      }else if(remainder(iteration,keep)==0){
        if(firstTime){
          oldbetasSum=oldBetas;
          firstTime=false; 
        }else{
          oldbetasSum+=oldBetas;
          oldbetasN+=1.0;      
        }
      }
      
      
    }
    
  }
  return Rcpp::List::create(
    Rcpp::Named("oldbetasSum") = oldbetasSum,
    Rcpp::Named("convergeStartedAt") = convergeStartedAt
  );
  
}