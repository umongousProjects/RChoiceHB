#include "HBfunctions.h"

void rGibbs(const mat &rGibbs_y,const  int nu,const  mat V,const  int nPar,const  int nRespondents,mat &rmultireg_B, 
                mat &rooti, Environment &myEnv,double IRd) {
        
        mat rGibbs_A = mat(1, 1);
        rGibbs_A(0, 0) = 0.01;
        mat rGibbs_X = mat(nPar, 1, fill::ones);
        //rmultireg
        mat rGibbs_RA = mat(1, 1);
        rGibbs_RA(0, 0) = 0.1;
        mat rGibbs_W = mat(nRespondents + 1, 1, fill::ones);
        rGibbs_W(nRespondents, 0) = 0.1;
        mat rGibbs_Z = rGibbs_RA * rmultireg_B;
        rGibbs_Z = join_cols(rGibbs_y, rGibbs_Z);
        mat IR = mat(1, 1);
        IR(0, 0) = IRd;
        mat Btilde = rGibbs_W.t() * rGibbs_Z;
        Btilde = (IR.t()*IR) * Btilde;
        
        mat rGibbs_S = rGibbs_W * Btilde;
        rGibbs_S = rGibbs_Z - rGibbs_S;
        rGibbs_S = rGibbs_S.t() * rGibbs_S;

        mat VpS = V + rGibbs_S;
        Function chol2inv2 = myEnv["chol2inv2"];
        Function rwishart = myEnv["rwishart"];
        Function backSolve_rGibbs = myEnv["backSolve_rGibbs"];
        
        mat cholVpS = as < mat > (chol2inv2(Rcpp::Named("x", VpS)));
                                
        //rwishart
        int rwishart_nu = nu + nRespondents;
        Rcpp::List res = rwishart(Rcpp::Named("nu", rwishart_nu), Rcpp::Named("V", cholVpS));
        rmultireg_B = Btilde;
        //arma::arma_rng::set_seed_random();
        mat randNums(1, nPar, fill::randn);
        randNums = IR * randNums;
        mat CI = as < mat > (res["CI"]);
        randNums *=  CI;
        rmultireg_B += randNums;
        mat rmultireg_IW = as < mat > (res["IW"]);
        //back to rGibbs
        rooti = as < mat > (backSolve_rGibbs(Rcpp::Named("sigm", rmultireg_IW)));
}

mat doV(const int nPar,const  int numOfVars,const  IntegerVector &featureLengths,const  int nu) {
        mat V = mat(nPar, nPar, fill::zeros); //not complete,see 135+
        int v_ind = 0;
        double v=1.0;
        double pcov = 0;
        double pvar = 0;
        mat temp;

        for (int i = 0; i < numOfVars; i++) {
                pcov = -1.0 / featureLengths(i);
                pvar = (featureLengths(i) - 1.0) / featureLengths(i);
                temp = eye < mat > (featureLengths(i) - 1, featureLengths(i) - 1);
                temp = pvar * temp;
                for (int ii = 0; ii < featureLengths(i) - 1; ii++) {
                        for (int iii = 0; iii < featureLengths(i) - 1; iii++) {
                                if (ii != iii) {
                                        temp(ii, iii) = pcov;
                                }
                        }
                }
                int ii_i = 0;
                int iii_i;
                for (int ii = v_ind; ii < v_ind + featureLengths(i) - 1; ii++) {
                        iii_i = 0;
                        for (int iii = v_ind; iii < v_ind + featureLengths(i) - 1; iii++) {
                                V(ii, iii) = temp(ii_i, iii_i);
                                iii_i += 1;
                        }
                        ii_i += 1;
                }
                v_ind = v_ind + featureLengths(i) - 1;
        }
        V = V * nu;
        V = V * v;
        return V;
}
mat doX(const int nPar,const  int numOfVars,const  IntegerVector &featureLengths,const  mat &X_cat) {
        mat X = mat(X_cat.n_rows, nPar, fill::zeros);
        int col_ind = 0;
        for (int thisVar = 0; thisVar < numOfVars; thisVar++) {
                for (int thisLevel = 1; thisLevel < featureLengths(thisVar); thisLevel++) {
                        for (int i = 0; i < X_cat.n_rows; i++) {
                                if (X_cat(i, thisVar) == thisLevel) {
                                        X(i, col_ind) = 1;
                                }
                                if (X_cat(i, thisVar) == featureLengths(thisVar)) {
                                        X(i, col_ind) = -1;
                                }
                        }
                        col_ind++;
                }
        }
        return X;
}

mat getLLMnl(const mat &X,const  int nPar,const  mat &newbeta,const  int packagesPerScreen,const  uvec &yIndices,
        const  int numberOfScreens) {
        mat repetitionColumn = mat(packagesPerScreen*numberOfScreens, 1, fill::ones); //
        mat tempXbeta = kron(newbeta, repetitionColumn);
        tempXbeta %= X;
        mat repetitionColumn2 = mat(nPar, 1, fill::ones);
        tempXbeta *= repetitionColumn2;

        mat Xbeta = tempXbeta;
        Xbeta.reshape(packagesPerScreen, Xbeta.n_rows / packagesPerScreen);

        mat xByY = tempXbeta.elem(yIndices);
        xByY.reshape(numberOfScreens, tempXbeta.n_rows / (numberOfScreens * packagesPerScreen));

        Xbeta.for_each([](mat::elem_type & val) {
                val = exp(val);
        });

        mat repetitionColumn3 = mat(1, packagesPerScreen, fill::ones);
        mat denom = repetitionColumn3 * Xbeta;
        denom.for_each([](mat::elem_type & val) {
                val = log(val);
        });
        denom.reshape(numberOfScreens, denom.n_cols / numberOfScreens);

        mat repetitionColumn4 = mat(1, numberOfScreens, fill::ones);
        mat ll = xByY - denom;
        ll = repetitionColumn4 * ll;

        return (ll.t());
}

mat getLndMvn(const mat &x,const mat &mu,const mat &rooti,const  int nPar){
        mat tempXminMu=x;
        tempXminMu.each_row() -=mu;
        mat z=tempXminMu*rooti;
        vec mainDiag = rooti.diag();
        double sumMainDiag=0;
        mat::iterator a = mainDiag.begin();
        mat::iterator b = mainDiag.end();
        for(mat::iterator i=a; i!=b; ++i)
        {
        sumMainDiag+=log(*i);
        }
        mat sumZSquared=z%z;
        mat repetitionColumn=mat(nPar,1,fill::ones);
        sumZSquared*=repetitionColumn;
        sumZSquared*=-0.5;
        double firstAndLastElement=-(nPar/2.0)*log(2.0*M_PI)+sumMainDiag;
        mat logs=sumZSquared+firstAndLastElement;
        return(logs);    
}

void mnlRwMetropOnce(const uvec &yIndices,const mat &X, mat &oldbeta,mat &oldll, 
        double s,const  mat &incroot, const  mat &betabar,const  mat &rootpi,const  int nRespondents,const  int numberOfScreens, 
        const int packagesPerScreen, const uvec &yIndices_old,const int nPar, int &stay) {

        //arma::arma_rng::set_seed_random();               
        mat randomNumbers = randn(nRespondents,nPar);
        mat increment=randomNumbers*incroot;
        increment*=s;
        mat newbeta=oldbeta+increment;
        mat newll=getLLMnl(X,nPar,newbeta,packagesPerScreen,yIndices_old,numberOfScreens);
        mat newlpost=newll+getLndMvn(newbeta,betabar,rootpi,nPar);

        mat ldiff=newlpost-oldll-getLndMvn(oldbeta,betabar,rootpi,nPar);
        mat alpha=ldiff;
        alpha.for_each( [](mat::elem_type& val) { val = exp(val); } );
        alpha.elem( find(alpha > 1.0) ).ones();
        //arma::arma_rng::set_seed_random();
        vec  unif = randu<vec>(nRespondents);
        unif.elem( find(alpha == 1.0) ).zeros();
        uvec good=find(unif <= alpha);

        unsigned int *a = good.begin();
        unsigned int *b = good.end();
        for(unsigned int * i=a; i!=b; ++i)
        {
                oldbeta.row(*i)=newbeta.row(*i);
                oldll.row(*i)=newll.row(*i);
        }
        stay=nRespondents- good.n_elem;
}