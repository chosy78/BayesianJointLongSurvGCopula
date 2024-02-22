//[[Rcpp::depends(RcppArmadillo)]]
#include <RcppArmadillo.h>
#include "Gcopula_jointMAR.h"

using namespace Rcpp;

//[[Rcpp::export]]
List MCMC_pwe(List & priors, arma::mat & Y, arma::mat & X, arma::mat & Z, arma::vec & v, arma::vec & obsind, arma::vec & Tst, arma::vec & delta, arma::mat & dat, arma::cube & W, List & inits, arma::vec & timeint, List & dim, List & slice_param, double & MAR, double & structured, List & IMidx, int & M, int & burnin, int & MHmcmc, int & thin, int & saving, double & accepttol, String & filename) {

  //##---------------------- PRIORS -------------------------##
Rcout<<"Priors..\n";  
  double a0, b0, c0, d0, sigr2, ar11, ar12, ar1priorBETA ; 
  arma::vec alpha0 ; arma::mat Sigalpha0 ;
  
  a0    = priors["a0"]    ; b0       = priors["b0"]       ;
  c0    = priors["c0"]    ; d0       = priors["d0"]       ;
  sigr2 = priors["sigr2"] ;
  arma::vec mu0   = priors["mu0"]   ; arma::mat Sig0     = priors["Sig0"]     ;
  arma::vec beta0 = priors["beta0"] ; arma::mat Sigbeta0 = priors["Sigbeta0"] ;
  if (structured==1.0) {
    ar11 = priors["ar11"]; ar12 = priors["ar12"]; ar1priorBETA = priors["ar1priorBETA"] ;
    
    arma::vec alpha0tmp    = priors["alpha0"]   ; alpha0    = alpha0tmp; 
    arma::mat Sigalpha0tmp = priors["Sigalpha0"]; Sigalpha0 = Sigalpha0tmp;
  }

  double consto, constmu, constalpha ; 
  consto = 1.0 ; constmu = 1.0 ; constalpha = 1.0 ; 
  //##-------------------------------------------------------##

  //##---------------------------- DATA ------------------------------##
Rcout<<"Data input..\n";
  arma::mat Yorigin, Xorigin; Yorigin = Y ; Xorigin = X;
  int m, n, J,J1, MM, q, mj, corrlength, muind; m = Y.n_cols ; n = Y.n_rows; J = dim["J"] ; J1 = J-1 ; MM = dim["M"] ; corrlength = m*(m-1)/2 ; mj = m/J ;
  double Xtmax, Xtmin;

  arma::mat Xm    = X(arma::span::all,arma::span(1,m)) ;
  arma::vec Xt    = X(arma::span::all,0) ;
  arma::mat Xobs  = X.rows(arma::find(v==1)) ; 
  arma::mat Xobsm = Xm.rows(arma::find(v==1)) ; 
  arma::vec Xobst = Xt.elem(arma::find(v==1)) ; 
  arma::mat Xcen  = X.rows(arma::find(v==0));
  arma::mat Xcenm = Xm.rows(arma::find(v==0));
  arma::vec Xcent = Xt.elem(arma::find(v==0));
  Xtmax = max(Xt) ; Xtmin = min(Xt);

  double Wexist;
  if (W.size()>1) {
    Wexist = 1.0;
    muind = W.n_cols;
  } else {
    Wexist = 0.0;
    muind = m;
  }
  
  q = muind / J ; 

  arma::vec nobscen = obsind.elem(arma::find(v==0));
  arma::vec nobsobs = obsind.elem(arma::find(v==1));
  
  // Intermittent Missing
Rcout<<"Missing..\n";
  arma::uvec uvec0(1) ; uvec0.fill(0);
  List IMcolnum_list = IMidx["colnum"]; List IMuvec_list = IMidx["idx"]; List IMcenidx_list=IMidx["cenidx"] ; List IMobsidx_list=IMidx["obsidx"];
  arma::uvec IMuvecfullyobserved = IMidx["fullyobservedidx"]; IMuvecfullyobserved = IMuvecfullyobserved-1;
  arma::uvec IMcenfullyobserved = IMidx["cenfullyobservedidx"] ; arma::uvec IMobsfullyobserved = IMidx["obsfullyobservedidx"]; IMcenfullyobserved = IMcenfullyobserved-1; IMobsfullyobserved=IMobsfullyobserved-1;
  
  arma::uvec cenIM = arma::find(nobscen==MM); arma::uvec obsIM = arma::find(nobsobs==MM); 
  int cenIMn = cenIM.n_elem; int obsIMn = obsIM.n_elem;  
  //##----------------------------------------------------------------##

  //##------------------------------------------- INITS ----------------------------------------------------##
Rcout<<"Initials..\n";  
  arma::vec sig2   = inits["sig2"]   ; int sig2n = sig2.n_elem;
  arma::vec mu     = inits["mu"]     ; 
  arma::vec lambda = inits["lambda"] ; 
  arma::vec beta   = inits["beta"]   ;
  arma::vec r12    = inits["r12"]    ; 
  arma::mat Gamma  = inits["Gamma"]  ; arma::mat Gammamm = Gamma(arma::span(1,m), arma::span(1,m)) ; arma::mat Gammainv = inv(Gamma) ; arma::mat Gammamminv = inv(Gammamm) ; double Gammastr_ar1=inits["Gammastr_ar1"];
  arma::mat Gammatmp = Gamma;

  arma::vec ar1vec(1), alpha; double ar1; //arma::mat Sigma; 
  ar1 = inits["ar1"] ;
  arma::vec alphatmp = inits["alpha"] ; alpha = alphatmp ; 
  int K, zm ;     K = lambda.n_elem ; zm = beta.n_elem ;
  //##------------------------------------------------------------------------------------------------------##

  //##------------Y and copula variable transformation-----------------##
Rcout<<"Transformed Data..\n";  
  arma::mat Yimp = Y;
  arma::uvec missind(Yimp.size()); missind.fill(n*m*2); arma::uvec missind1; 

  arma::vec sig2d(m) ; arma::mat SSIG(n,m); arma::mat Sig2(n,J); 
  arma::mat ymu(n,m);

  if (MAR==1.0) {  
    if (sig2n==1) {
      double ss = sig2(0);
      sig2d.fill(ss);
    } else {
      for (int i=0; i<J; i++) {
        arma::uvec indices(mj) ; int nn = 0 ;
        for (int ii=i; ii<m; ii=ii+J ) {
          indices(nn) = ii;
          nn +=1;
        }
        arma::vec ss(mj) ; ss.fill(sig2(i));
        sig2d(indices) = ss;
      }
    }
    arma::mat vv = diagmat(sig2d);

    if (Wexist==1.0) {
      arma::mat Wmu(m,n), rnd(m,n); arma::vec mm(m,arma::fill::zeros); 
      for (int i=0; i<n; i++)
        Wmu(arma::span::all,i) = W.slice(i) * mu;
      rnd  = arma::mvnrnd(mm,vv,n);
      Yimp = Wmu.t() + rnd.t();
    } else {
      Yimp = arma::mvnrnd(mu,vv,n);
    }

    Yimp.elem(find_finite(Yorigin)) = Yorigin.elem(find_finite(Yorigin));
    for (int i=0; i<n; i++) {
      int indl = J*(mj-obsind(i));

      arma::rowvec nanvec(indl) ; nanvec.fill(arma::datum::nan);
      if (indl>0)
        Yimp(i,arma::span(obsind(i)*J,m-1)) = nanvec;
    }

    arma::uvec yimpfinite = find_finite(Yimp) ; arma::uvec yoriinfi = find_nonfinite(Yorigin);
    int yimpn, yoriginn ; yimpn = yimpfinite.n_elem; yoriginn = yoriinfi.n_elem;
  
    int nn = 0;
    for (int i=0; i<yimpn; i++) {
      double yyimp = yimpfinite(i);
      for (int j=0; j<yoriginn; j++) {
        if (yyimp==yoriinfi(j)) {
          missind(nn) = yyimp;
          nn +=1;
        }
      }
    }
    missind = missind(arma::find(missind<(n*m*2))); missind1 = missind+1;
  }

  // Survival Data //
  arma::vec cumsum,survind ; cumsum = survind = arma::zeros(n);    
  //##----------------------------------------------------------------##

  //##------------------ BAGS ------------------------##
Rcout<<"Make bags..\n";  
  arma::mat SIG2(sig2n,M, arma::fill::zeros);
  arma::mat MU(muind,M, arma::fill::zeros);
  int rncol;
  if (structured==1.0) {
    if (Gammastr_ar1==4)
      rncol = 1;
    else 
      rncol = J*(J-1)/2;
  } else { 
    rncol = (m+1)*m/2;
  }
  
  arma::mat R12(rncol,M, arma::fill::zeros);
  arma::mat LAMBDA(K,M, arma::fill::zeros);
  arma::mat BETA(zm,M, arma::fill::zeros);
  arma::mat AR1(1,M, arma::fill::zeros);
  arma::mat ALPHA(m,M, arma::fill::zeros);

  arma::mat R12MH(rncol,MHmcmc, arma::fill::zeros);
  arma::mat MuMH(muind,MHmcmc, arma::fill::zeros);
  arma::mat AlphaMH(m,MHmcmc, arma::fill::zeros);

  int mhind, mhindmu, mhindalpha; mhind=mhindmu=mhindalpha=0;
  arma::mat MHvarmu, MHvaralpha, MHvar ; arma::mat VV(muind,muind); arma::mat VValpha(alpha.n_elem,alpha.n_elem) ; arma::mat VVr12(r12.n_elem,r12.n_elem) ;
  //##------------------------------------------------##

  //##------------------------------- COMPUTING TIMES --------------------------------##
  arma::vec timemu(M*thin + burnin + 2*MHmcmc, arma::fill::zeros); arma::vec timemuu((M*thin + burnin + 2*MHmcmc)*muind, arma::fill::zeros);
  arma::vec timesig2, timegamma, timebeta, timelambda, timesigb2, timeb, timear1, timealpha; timesig2 = timegamma = timebeta = timelambda = timesigb2 = timeb = timear1 = timealpha = timemu;
  arma::vec timemuint,timemupro; timemuint = timemupro = timemuu;
  //##--------------------------------------------------------------------------------##

  double wlong = slice_param["w"]; double slice_w_big = 20.0; double slice_w_small = 0.3; 
  double slice_m = slice_param["m"];

  int mcmcind=0; 
  int totalMCMC = M*thin + burnin + 2*MHmcmc+1; int burnmh1 = burnin + MHmcmc; int burnmh = burnin+2*MHmcmc; 

  double accpt,accptmu,accptalpha,accptind ; accpt=accptmu=accptalpha=accptind=0 ; 

  // Lists for log-likelihoods //
  List idx = List::create(Named("nobsobs")=nobsobs, _["nobscen"]=nobscen, _["uvec0"]=uvec0,_["IMcolnum"]=IMcolnum_list,_["IMuvec"]=IMuvec_list,_["IMuvecfullyobserved"]=IMuvecfullyobserved,_["IMcenidx"]=IMcenidx_list,_["IMobsidx"]=IMobsidx_list,_["IMcenfullyobserved"]=IMcenfullyobserved,_["IMobsfullyobserved"]=IMobsfullyobserved);

  List params_missing = List::create(Named("Yorigin")=Yimp, _["missind"]=missind, _["v"]=v, _["Xt"]=Xt, _["Xcent"]=Xcent, _["W"]=W,_["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Gamma"]=Gamma,_["mu"]=mu,_["sig2"]=sig2,_["idx"]=idx);
  List params_mu      = List::create(Named("Y")=Y,                _["v"]=v,                   _["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Xt"]=Xt,                                                                          _["mu0"]=mu0, _["Sig0"]=Sig0, _["sig2"]=sig2,                 _["Gamma"]=Gamma,_["W"]=W,_["idx"]=idx);
  List params_sig2    = List::create(Named("Y")=Y,                _["v"]=v,                   _["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Xt"]=Xt,                                                                          _["c0"]=c0, _["d0"]=d0, _["mu"]=mu,                           _["Gamma"]=Gamma,_["W"]=W,_["idx"]=idx);
  List params_lambda  = List::create(Named("Time")=Tst, _["Z"]=Z, _["v"]=v, _["delta"]=delta, _["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Xtmax"]=Xtmax,_["Xtmin"]=Xtmin, _["Xm"]=Xm,_["Xcenm"]=Xcenm,_["timeint"]=timeint, _["a0"]=a0, _["b0"]=b0, _["beta"]=beta,                       _["Gamma"]=Gamma,_["idx"]=idx);
  List params_beta    = List::create(Named("Time")=Tst, _["Z"]=Z, _["v"]=v, _["delta"]=delta, _["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Xtmax"]=Xtmax,_["Xtmin"]=Xtmin, _["Xm"]=Xm,_["Xcenm"]=Xcenm,_["timeint"]=timeint, _["beta0"]=beta0, _["Sigbeta0"]=Sigbeta0, _["lambda"]=lambda, _["Gamma"]=Gamma,_["idx"]=idx);
  List params_ar1     = List::create(Named("Y")=Y, _["v"]=v,_["n"]=n,_["M"]=MM,_["J"]=J, _["Xobs"]=Xobs,_["Xcen"]=Xcen,_["Xcenm"]=Xcenm,_["Xcent"]=Xcent, _["r12"]=r12, _["alpha"]=alpha, _["ar11"]=ar11,_["ar12"]=ar12,_["idx"]=idx,_["ar1priorBETA"]=ar1priorBETA,_["cenIM"]=cenIM,_["obsIM"]=obsIM,_["Gammastr_ar1"]=Gammastr_ar1);
  List params_r12uns  = List::create(Named("Y")=Y, _["n"]=n,_["m"]=m,_["M"]=MM,_["J"]=J, _["Xobs"]=Xobs,_["Xcen"]=Xcen,_["Xcenm"]=Xcenm,_["Xcent"]=Xcent, _["sigr2"]=sigr2, _["idx"]=idx,_["cenIM"]=cenIM,_["obsIM"]=obsIM);
  List params_alpha   = List::create(Named("Y")=Y, _["v"]=v,_["n"]=n,_["M"]=MM,_["J"]=J, _["Xobs"]=Xobs,_["Xcen"]=Xcen,_["Xcenm"]=Xcenm,_["Xcent"]=Xcent, _["r12"]=r12,_["ar1"]=ar1,_["alpha0"]=alpha0,_["Sigalpha0"]=Sigalpha0,_["mh"]=0.0,_["idx"]=idx,_["cenIM"]=cenIM,_["obsIM"]=obsIM,_["Gammastr_ar1"]=Gammastr_ar1);
  List params_r12str  = List::create(Named("Y")=Y, _["v"]=v,_["n"]=n,_["M"]=MM,_["J"]=J, _["Xobs"]=Xobs,_["Xcen"]=Xcen,_["Xcenm"]=Xcenm,_["Xcent"]=Xcent, _["sigr2"]=sigr2, _["ar1"]=ar1,_["alpha"]=alpha,_["idx"]=idx,_["cenIM"]=cenIM,_["obsIM"]=obsIM,_["Gammastr_ar1"]=Gammastr_ar1);

  String param_type_missing = "missing"; String param_type_mu = "mu"; String param_type_sig2 = "sig2"; String param_type_sig2mult = "sig2mult"; String param_type_lambda = "lambda"; String param_type_beta = "beta"; String param_type_ar1 = "ar1"; String param_type_r12uns = "Gamma"; String param_type_alpha = "alpha"; String param_type_r12str = "Gammastr";
  double lower_inf = -arma::datum::inf; double lower0 = 0.001 ; double lower_neg1 = -0.99 ; double upper = arma::datum::inf ; double upper1 = 0.99;
  // ------------------------- //
      
  //## ------------- ##
  //## MCMC SAMPLING ##
  //## ------------- ##
Rcout<<"Start MCMC sampling..\n";  
  for (int mcmc=0; mcmc<totalMCMC; mcmc++) {
    if (mcmc % 1000 == 0)
      Rcout << mcmc << "\n"; 

    // ------------ //
    // ---- MU ---- //
    // ------------ //
    if (MHmcmc == 0 || mcmc <= burnmh1) {      
      params_mu["Y"] = Y; params_mu["Xt"] = Xt; params_mu["sig2"] = sig2; params_mu["Gamma"] = Gamma; 
      mu = slice_sample_cpp(param_type_mu,mu,slice_w_big,slice_m,lower_inf,upper, params_mu);
    } else {
      for (int muvec=0; muvec<J; muvec++) {
        arma::vec mm = mu(arma::span(q*muvec, (q*(muvec+1))-1));
        arma::vec muparnew = arma::mvnrnd(mm,VV(arma::span(q*muvec, (q*(muvec+1))-1),arma::span(q*muvec, (q*(muvec+1))-1)),1);
        arma::vec munew    = mu;                    munew(arma::span(q*muvec, (q*(muvec+1))-1)) = muparnew;

        int tmpint = arma::datum::nan; double tmp = arma::datum::nan; double tmpmh = 1.0;
        double lpold = mupostdist_cpp(Y, v, n, m, MM, J, Xt, tmp, mu, tmpint, mu0, Sig0, sig2, Gamma, W, tmpmh, idx);
        double lpnew = mupostdist_cpp(Y, v, n, m, MM, J, Xt, tmp, munew, tmpint, mu0, Sig0, sig2, Gamma, W, tmpmh, idx);
        
        arma::vec tmp1 = arma::ones(2); tmp1(1) = exp(lpnew-lpold); 
        double accprob = arma::min(tmp1);

        double runif = arma::randu<double>();
        if (runif < accprob) {
          accptmu = accptmu + 1/J;
          mu = munew;
        }
      }  
     
      if (mcmc == burnmh) {
        double tmp;
        if (accptmu/burnmh < 0.01)
          tmp = 1e-30;
        else 
          tmp = accptmu/burnmh;
        constmu = (constmu*R::qnorm(accepttol/2, 0.0,1.0,TRUE,FALSE)) / R::qnorm(0.5*tmp, 0.0,1.0,TRUE,FALSE);

        Rcout << "accprob:" << accptmu/burnmh << "\n" ;
        Rcout << "const:" << constmu << "\n" ;

        VV = (constmu*constmu) * MHvarmu;
      }
    }
   
    if (MHmcmc > 0) {
      if ((mcmc > burnin) && (mcmc <= burnmh1)) {
        MuMH(arma::span::all,mhindmu) = mu;
        mhindmu +=1;
      }
      if (mcmc==burnmh1) {
        MHvarmu  = arma::cov(MuMH.t());
        Rcout << "MHvarmu" << MHvarmu << "\n";
      }
    }

    // -------------- //
    // ---- SIG2 ---- //
    // -------------- //
    params_sig2["Y"] = Y; params_sig2["Xt"] = Xt; params_sig2["mu"] = mu; params_sig2["Gamma"] = Gamma; 
    if (sig2n==1) 
      sig2 = slice_sample_cpp(param_type_sig2,sig2,slice_w_big,slice_m,lower0,upper, params_sig2);
    else 
      sig2 = slice_sample_cpp(param_type_sig2mult,sig2,slice_w_big,slice_m,lower0,upper, params_sig2);
    
    if (sig2n==1) {
      double ss = sig2(0);
      
      SSIG.fill(ss);
    } else {
      Sig2.each_row() = reshape(sig2,1,J);

      for (int i=0;i<MM;i++)
        SSIG(arma::span::all,arma::span(J*i,J*(i+1)-1)) = Sig2;  
    }

    if (Wexist==1.0) {
      arma::mat Wmu(m,n);
      for (int i=0; i<n; i++)
        Wmu(arma::span::all,i) = W.slice(i) * mu;
      ymu = Y - Wmu.t();
    } else {
      //arma::mat Yt = Y.t();
      arma::mat ymut(m,n);
      for (int i=0; i<n; i++)
        ymut(arma::span::all,i) = Y(i,arma::span::all).t() - mu;
      ymu = ymut.t();
    }

    Xm = ymu / sqrt(SSIG);
    for (int i=0; i<m; i++)
      Xm(arma::span::all,i) = Xm(arma::span::all,i) * sqrt(Gamma(i+1,i+1));

    Xobsm = Xm.rows(arma::find(v==1));
    Xcenm = Xm.rows(arma::find(v==0));

    Xtmax = max(Xt); Xtmin = min(Xt);
    
    // ---------------- //
    // ---- LAMBDA ---- //
    // ---------------- //
    params_lambda["Xtmax"]=Xtmax; params_lambda["Xtmin"]=Xtmin; params_lambda["Xm"] = Xm; params_lambda["Xcenm"] = Xcenm; params_lambda["beta"] = beta; params_lambda["Gamma"] = Gamma; 
    lambda = slice_sample_cpp(param_type_lambda,lambda,wlong,slice_m,lower0,upper, params_lambda);
    
    // -------------- //
    // ---- BETA ---- //
    // -------------- //
    params_beta["Xtmax"]=Xtmax; params_beta["Xtmin"]=Xtmin; params_beta["Xm"] = Xm; params_beta["Xcenm"] = Xcenm; params_beta["lambda"] = lambda; params_beta["Gamma"] = Gamma; 
    beta = slice_sample_cpp(param_type_beta,beta,wlong,slice_m,lower_inf,upper, params_beta);
    
    // ----------------------- //
    // ---- Survival Data ---- //
    // ----------------------- //
    for (int i=0; i<n; i++) {
      if (delta(i)==1) 
        cumsum(i) = 0;
      else
        cumsum(i) = sum(lambda(arma::span(0,delta(i)-2)) % (timeint(arma::span(1,delta(i)-1)) - timeint(arma::span(0,delta(i)-2))));
      survind(i) = exp(-( lambda(delta(i)-1) * (Tst(i)-timeint(delta(i)-1)) + cumsum(i) ) * exp(accu(Z(i,arma::span::all) % beta.t())));
    }

    for (int i=0; i<n; i++) {
      Xt(i) = R::qnorm(1-survind(i), 0.0, 1.0, TRUE, FALSE);
      if (!arma::is_finite(Xt(i))) {  
        if (survind(i)<0.1)
          Xt(i) = Xtmax;
        else 
          Xt(i) = Xtmin;
      }
    }
    X(arma::span::all,0) = Xt;
    X(arma::span::all,arma::span(1,m)) = Xm;
    Xobs  = X.rows(arma::find(v==1)) ; 
    Xobst = Xt.elem(arma::find(v==1)) ; 
    Xcen  = X.rows(arma::find(v==0));
    Xcent = Xt.elem(arma::find(v==0));

    // --------------------- //
    // ---- AR1 & ALPHA ---- //
    // --------------------- //
    if (structured==1.0) {
      // ------- //
      // - AR1 - //
      // ------- //
      
      params_ar1["Y"] = Y; params_ar1["Xobs"] = Xobs; params_ar1["Xcen"] = Xcen; params_ar1["Xcenm"] = Xcenm; params_ar1["Xcent"] = Xcent; params_ar1["r12"] = r12; params_ar1["alpha"] = alpha; 
      ar1vec.fill(ar1);
      if (ar1priorBETA==1.0) {
        ar1vec = slice_sample_cpp(param_type_ar1,ar1vec,slice_w_small,slice_m,lower0,upper1, params_ar1);
        ar1 = ar1vec(0);
      } else {
        ar1vec = slice_sample_cpp(param_type_ar1,ar1vec,slice_w_small,slice_m,lower_neg1,upper1, params_ar1);
        ar1 = ar1vec(0);
      } 
    }

    // --------- //
    // - GAMMA - //
    // --------- //
    if (MHmcmc ==0 || mcmc <= burnmh1) {

      if (structured==0.0) {
        params_r12uns["Y"] = Y; params_r12uns["Xobs"] = Xobs; params_r12uns["Xcen"] = Xcen; params_r12uns["Xcenm"] = Xcenm; params_r12uns["Xcent"] = Xcent; 
        r12 = slice_sample_cpp(param_type_r12uns,r12,wlong,slice_m,lower_inf,upper, params_r12uns);

        Gammatmp = Gammacal(r12,m) ; Gamma = Gammatmp ; 
        for (int ii=0; ii<m+1; ii++) {
          Gamma(arma::span::all,ii) /= sqrt(Gammatmp(ii,ii)) ;
          Gamma(ii,arma::span::all) /= sqrt(Gammatmp(ii,ii)) ;
        }
        Gammamm = Gamma(arma::span(1,m), arma::span(1,m)) ; Gammainv = inv(Gamma) ; Gammamminv = inv(Gammamm) ;
 
      } else {
        // ---------------- //
        // - ALPHA & Gamma- //
        // ---------------- //
        params_alpha["Y"] = Y; params_alpha["Xobs"] = Xobs; params_alpha["Xcen"] = Xcen; params_alpha["Xcenm"] = Xcenm; params_alpha["Xcent"] = Xcent; params_alpha["r12"] = r12; params_alpha["ar1"] = ar1; 
        alpha = slice_sample_cpp(param_type_alpha,alpha,wlong,slice_m,lower_inf,upper, params_alpha);
        if (J>1) {
          arma::mat Sigma; 
          if (Gammastr_ar1==3) {
            Sigma = AR1cal(ar1,m,Gammastr_ar1);
          } else if (Gammastr_ar1==4) {
            params_ar1["Y"] = Y; params_ar1["Xobs"] = Xobs; params_ar1["Xcen"] = Xcen; params_ar1["Xcenm"] = Xcenm; params_ar1["Xcent"] = Xcent; params_ar1["r12"] = ar1vec; params_ar1["alpha"] = alpha; params_ar1["ar1priorBETA"]=0.0;
            r12 = slice_sample_cpp(param_type_ar1,r12,slice_w_small,slice_m,lower_neg1,upper1, params_ar1);
            Sigma = kron(AR1cal(ar1,MM,Gammastr_ar1),AR1cal(r12(0),J,Gammastr_ar1)); 
          } else {
            params_r12str["Y"] = Y; params_r12str["Xobs"] = Xobs; params_r12str["Xcen"] = Xcen; params_r12str["Xcenm"] = Xcenm; params_r12str["Xcent"] = Xcent; params_r12str["alpha"] = alpha; 
            r12 = slice_sample_cpp(param_type_r12str,r12,wlong,slice_m,lower_inf,upper, params_r12str);
            Sigma = kron(AR1cal(ar1,MM,Gammastr_ar1),Gammacal(r12,J1)); 
          }
          Gamma(0,0) = 1.0; Gamma(0,arma::span(1,m)) = alpha.t(); Gamma(arma::span(1,m),0) = alpha; Gamma(arma::span(1,m),arma::span(1,m)) = Sigma; 
          Gammatmp = Gamma;  
         } 
        Gammamm = Gamma(arma::span(1,m), arma::span(1,m)) ; Gammainv = inv(Gamma) ; Gammamminv = inv(Gammamm) ;
      }
    } else {

      if (structured==1.0) {
        arma::vec alphanew = arma::mvnrnd(alpha,VValpha,1);

        int tmpint = arma::datum::nan; double tmp = arma::datum::nan; double tmpmh = 1.0;
        double oldalpha = alphapostdist_cpp(Y, v, n,MM,J, Xobs,Xcen,Xcenm,Xcent, r12,ar1,tmp, alpha, tmpint, alpha0,Sigalpha0, tmpmh, idx,cenIM,obsIM,Gammastr_ar1);
        double newalpha = alphapostdist_cpp(Y, v, n,MM,J, Xobs,Xcen,Xcenm,Xcent, r12,ar1,tmp, alphanew, tmpint, alpha0,Sigalpha0, tmpmh, idx,cenIM,obsIM,Gammastr_ar1);
        
        arma::vec tmp1 = arma::ones(2); tmp1(1) = exp(newalpha-oldalpha); 
        double accprob = arma::min(tmp1);

        double runif = arma::randu<double>();
        if (runif < accprob) {
          accptalpha +=1;
          alpha = alphanew;
        }
      
        if (mcmc == burnmh) {
          double tmp;
          if (accptalpha/burnmh < 0.01)
            tmp = 1e-30;
          else 
            tmp = accptalpha/burnmh;
          constalpha = (constalpha*R::qnorm(accepttol/2, 0.0,1.0,TRUE,FALSE)) / R::qnorm(0.5*tmp, 0.0,1.0,TRUE,FALSE);
          Rcout << "Alpha accprob:" << accptalpha/burnmh << "\n" ;
          Rcout << "Alpha const:" << constalpha << "\n" ;

          VValpha = constalpha*constalpha * MHvaralpha;
        }
      }

      arma::vec r12new = r12;
      if (!(structured==1.0 && J==1)) {
        r12new = arma::mvnrnd(r12,VVr12,1);
      }

      int tmpint = arma::datum::nan; double tmp = arma::datum::nan; double tmpmh = 1.0;
      List old, neww;
      if (structured==0.0) {        
        old = r12postdist_cpp(Y,n,m,MM,J, Xobs,Xcen,Xcenm,Xcent,tmp, r12, tmpint, sigr2, tmpmh, idx,cenIM,obsIM);
        neww = r12postdist_cpp(Y,n,m,MM,J, Xobs,Xcen,Xcenm,Xcent,tmp, r12new, tmpint, sigr2, tmpmh, idx,cenIM,obsIM);
      } else if (structured==1.0 && J!=1) {
        old = Gammapostdist_cpp(Y,v,n,MM,J, Xobs,Xcen,Xcenm,Xcent,tmp, r12, tmpint, sigr2, ar1, alpha, tmpmh, idx,cenIM,obsIM,Gammastr_ar1);
        neww = Gammapostdist_cpp(Y,v,n,MM,J, Xobs,Xcen,Xcenm,Xcent,tmp, r12new, tmpint, sigr2, ar1, alpha, tmpmh, idx,cenIM,obsIM,Gammastr_ar1);
      } 

      if (!(structured==1.0 && J==1)) { 
      
        double lpold = old["logpost"] ; arma::mat Gammatmp = old["Gamma"] ; Gamma = Gammatmp;
        for (int ii=0; ii<m+1; ii++) {
          Gamma(arma::span::all,ii) /= sqrt(Gammatmp(ii,ii)) ;
          Gamma(ii,arma::span::all) /= sqrt(Gammatmp(ii,ii)) ;
        }
        
        Gammamm = Gamma(arma::span(1,m), arma::span(1,m)) ; Gammainv = inv(Gamma) ; Gammamminv = inv(Gammamm) ;
        double lpnew = neww["logpost"] ; arma::mat Gammanew = neww["Gamma"] ; Gammatmp = Gammanew;
        for (int ii=0; ii<m+1; ii++) {
          Gammanew(arma::span::all,ii) /= sqrt(Gammatmp(ii,ii)) ;
          Gammanew(ii,arma::span::all) /= sqrt(Gammatmp(ii,ii)) ;
        }
        
        arma::mat Gammammnew = Gammanew(arma::span(1,m), arma::span(1,m)) ; arma::mat Gammainvnew = inv(Gammanew) ; arma::mat Gammamminvnew = inv(Gammammnew) ;
  
        arma::vec tmp1 = arma::ones(2); tmp1(1) = exp(lpnew-lpold); 
        double accprob = arma::min(tmp1);

        double runif = arma::randu<double>();
        if (runif < accprob) {
          accpt +=1;
          r12 = r12new;
          Gamma = Gammanew ; Gammamm = Gammammnew ; Gammainv = Gammainvnew ; Gammamminv = Gammamminvnew ;
        }
  
        if (mcmc == burnmh) {
          double tmp;
          if (accpt/burnmh < 0.01)
            tmp = 1e-30;
          else 
            tmp = accpt/burnmh;
          consto = (consto*R::qnorm(accepttol/2, 0.0,1.0,TRUE,FALSE)) / R::qnorm(0.5*tmp, 0.0,1.0,TRUE,FALSE);
          Rcout << "Gamma accprob:" << accpt/burnmh << "\n" ;
          Rcout << "Gamma const:" << consto << "\n" ;

          VVr12 = consto*consto*MHvar;
        }
      } 
    }

    if (J==1) {
      r12 = arma::zeros(1) ; int MM1 = MM+1;
      Gamma = AR1cal(ar1, MM1,Gammastr_ar1);
      Gammamm = Gamma(arma::span(1,m), arma::span(1,m)) ; Gammainv = inv(Gamma) ; Gammamminv = inv(Gammamm) ;
    }

    if (MHmcmc > 0) {
      if ((mcmc > burnin) && (mcmc <= burnmh1)) {
        R12MH(arma::span::all,mhind) = r12;
        AlphaMH(arma::span::all,mhind) = alpha;
        mhind +=1;
      }
      if (mcmc==burnmh1) {
        MHvar  = arma::cov(R12MH.t()); 
        MHvaralpha = arma::cov(AlphaMH.t());
        Rcout << "MHvar" << MHvar << "\n";
        Rcout << "MHvaralpha" << MHvaralpha << "\n";
      }  
    }

    int mcmcburnmh = mcmc-burnin-2*MHmcmc;
    if ((mcmc > burnmh) && (mcmcburnmh % thin == 0)) {
      MU(arma::span::all,mcmcind) = mu;
      SIG2(arma::span::all,mcmcind) = sig2;
      LAMBDA(arma::span::all,mcmcind) = lambda;
      BETA(arma::span::all,mcmcind) = beta;
      R12(arma::span::all,mcmcind) = r12;
      if (structured==1.0) {
        AR1(arma::span::all,mcmcind) = ar1;
        ALPHA(arma::span::all,mcmcind) = alpha;
      }
      mcmcind +=1;
    }
  }
  List Result = List::create(Named("M")=M, _["MU"]=MU, _["SIG2"]=SIG2, _["LAMBDA"]=LAMBDA, _["BETA"]=BETA, _["R12"]=R12, _["AR1"]=AR1, _["ALPHA"]=ALPHA, _["dat"]=dat,_["W"]=W, _["accpt"]=List::create(Named("accpt_gamma")=accpt/(M*thin + MHmcmc), _["accpt_mu"] = accptmu/(M*thin + MHmcmc), _["accpt_alpha"]=accptalpha/(M*thin + MHmcmc)));   
  return Result;

}
