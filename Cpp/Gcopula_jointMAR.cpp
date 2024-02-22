# define ARMA_DONT_USE_WRAPPER
//[[Rcpp::depends(RcppArmadillo)]]
//[[Rcpp::depends(RcppEigen)]]
//[[Rcpp::depends(RcppNumerical)]]
#include <RcppArmadillo.h>
#include <RcppNumerical.h>   
#include <vector>

using namespace Rcpp;
using namespace Numer;

/******************/
/* BASE FUNCTIONS */
/******************/
//[[Rcpp::export]]
arma::mat Gammacal(arma::vec & r12, int & m) {
  
  arma::mat R = arma::mat(m+1,m+1, arma::fill::eye);
  
  int ind = 0;
  for (int j = 0; j < m+1; j++)
  {
    for (int i = 0; i < m+1; i++) 
    {
      if (i < j)
      {
        R(i,j) = r12[ind];
        ind ++;
      }
    }	
  }
  
  arma::mat Sigma = inv(R.t() * R);
  arma::vec diag  = 1/sqrt(Sigma.diag());
  arma::mat Diag  = diagmat(diag);
  arma::mat Gamma = Diag * Sigma * Diag;
  
  return Gamma;
}

//[[Rcpp::export]]
arma::mat AR1cal(double & ar1, int & m, double & Gammastr_ar1) {
  arma::mat R = arma::mat(m,m,arma::fill::zeros);
  arma::mat RR = arma::mat(m,m,arma::fill::zeros);
  arma::mat Rho = arma::mat(m,m) ; Rho.fill(ar1);
  arma::mat exponent;
  
  if (Gammastr_ar1==1.0) {
    for (int i=1; i<m; i++) {
      arma::vec a(m) ; arma::rowvec b(m) ; a.fill(i); b.fill(i);
      R.col(i) = a;
      RR.row(i) = b;
    }
    exponent = exp(abs(R-RR)%log(Rho));
  } else {
    exponent = Rho; exponent.diag().ones();
  }
  return(exponent);
}


/********************************/
/* FIRST PAPER RELATED FUNCTION */
/********************************/
//[[Rcpp::export]]
double survcdflog_cpp(arma::mat & Xcenm, arma::vec & Xcent, int & m, int & n, int & M, int & J, arma::mat & Gamma, List & idx){ //arma::vec & nobs, arma::uvec & uvec0,arma::uvec & IMcolnum,arma::uvec & IMcenidx,arma::uvec & IMcenfullyobserved) {
	
  arma::vec nobs = idx["nobscen"] ; arma::uvec uvec0=idx["uvec0"] ; arma::uvec IMcenfullyobserved=idx["IMcenfullyobserved"];
  List IMcolnum_list=idx["IMcolnum"] ; List IMcenidx_list=idx["IMcenidx"];

	double result=0.0; 
	for (int ii=1; ii < M+1 ; ii++)
	{	  
		arma::mat Gammamt(ii*J,1);
		arma::mat Gammamminv(ii*J,ii*J);
		Gammamt    = Gamma(arma::span(1,ii*J),0);
		Gammamminv = inv(Gamma(arma::span(1,ii*J),arma::span(1,ii*J)));
		arma::mat Xm; arma::vec Xt; 
		Xm = Xcenm.rows(arma::find(nobs==ii)); Xt = Xcent(arma::find(nobs==ii));
    
		int nn = Xt.n_elem;
		if (nn > 0)
		{
      if (ii==M) {
        // FULLY OBSERVED
        arma::mat condmean = (Gammamt.t()*Gammamminv) * Xcenm.rows(IMcenfullyobserved).t();
        arma::mat condvar  = 1 - Gammamt.t()*(Gammamminv * Gammamt);
        arma::mat varvec   = arma::repelem(sqrt(condvar),1,IMcenfullyobserved.n_elem);

        arma::vec resulttmp(IMcenfullyobserved.n_elem); 
        for (int p=0; p<IMcenfullyobserved.n_elem; p++)
        {
          resulttmp(p) = 1-R::pnorm(Xcent(IMcenfullyobserved(p)),condmean(p), varvec(p), TRUE, FALSE);
          if (resulttmp(p)==0) resulttmp(p) = 1e-7;
        }
        result += sum(log(resulttmp));    

        // INTERMITTENT MISSING
        for (int i=0; i < IMcolnum_list.size(); i++){
          arma::uvec IMcolnum  = as<arma::uvec>(IMcolnum_list[i]); arma::uvec IMcenidx = as<arma::uvec>(IMcenidx_list[i]); IMcenidx = IMcenidx-1;
          arma::mat Gammamt(IMcolnum.n_elem,1); Gammamt = Gamma.submat(IMcolnum,uvec0);
          arma::mat Gammamminv(IMcolnum.n_elem,IMcolnum.n_elem); Gammamminv = inv(Gamma.submat(IMcolnum,IMcolnum));
          condmean = (Gammamt.t()*Gammamminv) * Xcenm.submat(IMcenidx,IMcolnum-1).t();
          condvar  = 1 - Gammamt.t()*(Gammamminv * Gammamt);
          varvec   = arma::repelem(sqrt(condvar),1,IMcenidx.n_elem);
          arma::vec resultttmp(IMcenidx.n_elem); 

          for (int p=0; p<IMcenidx.n_elem; p++)
          {
            resultttmp(p) = 1-R::pnorm(Xcent(IMcenidx(p)),condmean(p), varvec(p), TRUE, FALSE);
            if (resultttmp(p)==0) resultttmp(p) = 1e-7;
          }
          result += sum(log(resultttmp));    
        }
      } else {
        arma::mat condmean = (Gammamt.t()*Gammamminv) * Xm(arma::span::all,arma::span(0,ii*J-1)).t();
  			arma::mat condvar  = 1 - Gammamt.t()*(Gammamminv * Gammamt);
  			arma::mat varvec   = arma::repelem(sqrt(condvar),1,nn);
  			arma::vec resulttmp(nn); 
  			for (int p=0; p<nn; p++)
  			{
  				resulttmp(p) = 1-R::pnorm(Xt(p),condmean(p), varvec(p), TRUE, FALSE);
  				if (resulttmp(p)==0) resulttmp(p) = 1e-7;
  			}
  			result += sum(log(resulttmp));  
      }
	  }
  }
	return result;
}

//[[Rcpp::export]]
double obsliklog_cpp(arma::mat & Xobs, int & m, int & M, int & J, arma::mat & Gamma, List & idx){ //arma::vec & nobs, arma::uvec & IMcolnum0,arma::uvec & IMobsidx,arma::uvec & IMobsfullyobserved) {

  arma::vec nobs = idx["nobsobs"] ; arma::uvec IMobsfullyobserved=idx["IMobsfullyobserved"]; arma::uvec uvec0=idx["uvec0"]; 
  List IMcolnum_list=idx["IMcolnum"]; List IMobsidx_list=idx["IMobsidx"];

	double result=0.0; 
	for (int ii=1; ii < M+1 ; ii++)
	{
    arma::mat Xm; arma::mat Xobsx;
    Xm = Xobs.rows(arma::find(nobs==ii));
    
    if (ii==M) {
      // FULLY OBSERVED
      result += -0.5 * accu((Xobs.rows(IMobsfullyobserved) * (inv(Gamma)-arma::mat(m+1,m+1,arma::fill::eye))) % Xobs.rows(IMobsfullyobserved));

      // INTERMITTENT MISSING
      for (int i=0; i < IMcolnum_list.size(); i++){
        arma::uvec IMcolnum  = as<arma::uvec>(IMcolnum_list[i]); arma::uvec IMobsidx = as<arma::uvec>(IMobsidx_list[i]); IMobsidx = IMobsidx-1;
        arma::uvec IMcolnum0 = arma::join_cols(uvec0,IMcolnum);
        arma::mat Cov = inv(Gamma.submat(IMcolnum0,IMcolnum0)) - arma::mat(IMcolnum0.n_elem,IMcolnum0.n_elem,arma::fill::eye);
        Xobsx = Xobs.submat(IMobsidx,IMcolnum0);
        result += -0.5 * accu((Xobsx * Cov) % Xobsx);
      }

    } else {
  		arma::mat Gammamtinv = Gamma(arma::span(0,ii*J),arma::span(0,ii*J));
  		arma::mat Cov 		   = inv(Gammamtinv) - arma::mat(ii*J+1,ii*J+1, arma::fill::eye);
  		Xobsx   = Xm(arma::span::all, arma::span(0,ii*J));
  		result += -0.5 * accu((Xobsx * Cov) % Xobsx);
    }
	}
	return result;
}

//[[Rcpp::export]]
double censliklog_cpp(arma::mat & Xcenm, int & m, int & M, int & J, arma::mat & Gamma, List & idx) {

  arma::vec nobs = idx["nobscen"] ; arma::uvec IMcenfullyobserved=idx["IMcenfullyobserved"];
  List IMcolnum_list=idx["IMcolnum"] ; List IMcenidx_list=idx["IMcenidx"];

	double result=0.0; 
	for (int ii=1; ii < M+1 ; ii++)
	{
    arma::mat Xm; 
    arma::mat Xcenmx;
    Xm = Xcenm.rows(arma::find(nobs==ii));
    
    if (ii==M) {
      // FULLY OBSERVED
      result += -0.5 * accu((Xcenm.rows(IMcenfullyobserved) * (inv(Gamma(arma::span(1,M*J),arma::span(1,M*J)))-arma::mat(M*J,M*J,arma::fill::eye))) % Xcenm.rows(IMcenfullyobserved));  
      // INTERMITTENT MISSING
      for (int i=0; i < IMcolnum_list.size(); i++){
        arma::uvec IMcolnum  = as<arma::uvec>(IMcolnum_list[i]); arma::uvec IMcenidx = as<arma::uvec>(IMcenidx_list[i]); IMcenidx = IMcenidx-1;
        arma::mat Cov = inv(Gamma.submat(IMcolnum,IMcolnum)) - arma::mat(IMcolnum.n_elem,IMcolnum.n_elem,arma::fill::eye);
        Xcenmx = Xcenm.submat(IMcenidx,IMcolnum-1);
        result += -0.5 * accu((Xcenmx * Cov) % Xcenmx);
      }
    } else {
      arma::mat Cov = inv(Gamma(arma::span(1,ii*J),arma::span(1,ii*J))) - arma::mat(ii*J,ii*J, arma::fill::eye);
      Xcenmx  = Xm(arma::span::all, arma::span(0,ii*J-1));
      result += -0.5 * accu((Xcenmx * Cov) % Xcenmx);  
    }
	}
	return result;
}

//[[Rcpp::export]]
double missingliklog_cpp(arma::mat & Yorigin, double & ymiss, int & ind, arma::vec & missind, arma::vec & Ximp, arma::vec & v, arma::vec & Xt, arma::vec & Xcent, arma::cube & W, int & n, int & m, int & M, int & J, arma::mat & Gamma, arma::vec & mu, arma::vec & sig2, List & idx) {
  // MAR SHOULD NOT BE CONSIDERED, THUS NOT USED
  arma::mat Y = Yorigin;
  Y(missind(ind-1)) = ymiss;
  arma::mat ymu(n,m);
  if (W.size()>1) {
    arma::mat Wmu(m,n);
    for (int i=0; i<n; i++)
      Wmu(arma::span::all,i) = W.slice(i) * mu;
    ymu = Y - Wmu.t();
  } else {
    arma::mat Yt = Y.t();
    arma::mat ymut(m,n);
    for (int i=0; i<n; i++)
      ymut(arma::span::all,i) = Yt(arma::span::all,i) - mu;
    ymu = ymut.t();
  }

  arma::mat SSIG(n,m); 

  if (sig2.size()==1) {
    double ss = sig2(0);
    SSIG.fill(ss);
  } else {
    arma::mat Sig2(n,J); Sig2.each_row() = reshape(sig2,1,J);
    for (int i=0;i<M;i++)
      SSIG(arma::span::all,arma::span(J*i,J*(i+1)-1)) = Sig2;  
  }

  arma::mat Xmb = ymu / sqrt(SSIG); arma::mat Xm(n,m);
  arma::mat Xm2 = Xmb % Xmb;

  for (int i=0; i<m; i++)
    Xm(arma::span::all,i) = Xmb(arma::span::all,i) * sqrt(Gamma(i+1,i+1));

  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::mat Xcenm = Xm.rows(arma::find(v==0));

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);
  double datalik  = -0.5 * accu(Xm2.elem(find_finite(Xm2)));
  double result = cond + obs + cen + datalik;
  return result;
}
  

//[[Rcpp::export]]
Rcpp::List r12postdist_cpp(arma::mat & Y, int & n, int & m, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, double & r, arma::vec & r_old, int & ind, double & sigr2, double & mh, List & idx, arma::uvec & cenIM, arma::uvec & obsIM ) {

  arma::vec r12 = r_old;
  if (mh==0.0)
  	r12(ind-1) = r;

  arma::mat Gamma = Gammacal(r12, m);
  arma::vec nobsobs = idx["nobsobs"] ; arma::vec nobscen = idx["nobscen"]; 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double prior;
  if (mh==0.0)
  	prior = 0.5 * (1/sigr2) * r*r;
  else 
  	prior = 0.5 * (1/sigr2) * accu(r12 % r12);
  
  double gamobs = 0.0;
  double gamcen = 0.0;

  for (int i=1; i<M+1; i++)
  { 
    if (i==M) {
      for (int p=0; p<obsIM.n_elem; p++) { 
        arma::uvec ttmp = find_finite(Xobs(obsIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamobs += -0.5 * val;
      }
    } else {
      double val; double sign;  
      log_det(val, sign, Gamma(arma::span(0,i*J),arma::span(0,i*J)));
      gamobs += -0.5 * sum(nobsobs(arma::find(nobsobs==i)))/i * val;  
    }
  }

  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<cenIM.n_elem; p++) { 
        arma::uvec ttmp = 1+find_finite(Xcenm(cenIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamcen += -0.5 * val;
      }
    } else {
    	double val; double sign;
      log_det(val, sign, Gamma(arma::span(1,i*J),arma::span(1,i*J)));
    	gamcen += -0.5 * sum(nobscen(arma::find(nobscen==i)))/i * val;
    }
  }

  double logpost = gamobs + gamcen + cond + obs + cen - prior;
  
  List Result = List::create(Named("Gamma")=Gamma, _["logpost"]=logpost);

  return Result;
}

//[[Rcpp::export]]
Rcpp::List Gammapostdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, double & r, arma::vec & r_old, int & ind, double & sigr2, double & ar1, arma::vec & alpha, double & mh, List & idx, arma::uvec & cenIM, arma::uvec & obsIM, double & Gammastr_ar1 ) {

  int m = M*J;

  arma::vec r12 = r_old;
  if (mh==0.0)
    r12(ind-1) = r;

  int mm = M+1 ; int jj = J-1;

  arma::mat Sigma; 
  if (Gammastr_ar1==3) {
    Sigma = AR1cal(ar1,m,Gammastr_ar1);
  }
  else if (Gammastr_ar1==4) { // does not happen
    double r12double =r12(0);
    Sigma = kron(AR1cal(ar1,M,Gammastr_ar1),AR1cal(r12double,J,Gammastr_ar1)); 
  } else {
    Sigma = kron(AR1cal(ar1,M,Gammastr_ar1),Gammacal(r12,jj)); 
  }
  arma::mat Sigmainv = inv(Sigma);
  arma::mat Gammainv(m+1,m+1);
  Gammainv(0,0) = as_scalar(alpha.t() * Sigmainv * alpha) + 1;
  Gammainv(arma::span(1,m),arma::span(1,m)) = Sigmainv;
  Gammainv(0,arma::span(1,m)) = - alpha.t() * Sigmainv ; Gammainv(arma::span(1,m),0) = - Sigmainv * alpha;
  
  arma::mat Gamma = inv(Gammainv);
  arma::vec nobsobs = idx["nobsobs"] ; arma::vec nobscen = idx["nobscen"]; 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double prior;
  if (mh==0.0)
    prior = 0.5 * (1/sigr2) * r*r;
  else 
    prior = 0.5 * (1/sigr2) * accu(r12 % r12);

  double gamobs = 0.0; double gamcen = 0.0;
  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<obsIM.n_elem; p++) { 
        arma::uvec ttmp = find_finite(Xobs(obsIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamobs += -0.5 * val;
      }
    } else {
      double val; double sign;  
      log_det(val, sign, Gamma(arma::span(0,i*J),arma::span(0,i*J)));
      gamobs += -0.5 * sum(nobsobs(arma::find(nobsobs==i)))/i * val;  
    }
  }

  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<cenIM.n_elem; p++) { 
        arma::uvec ttmp = 1+find_finite(Xcenm(cenIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamcen += -0.5 * val;
      }
    } else {
      double val; double sign;
      log_det(val, sign, Gamma(arma::span(1,i*J),arma::span(1,i*J)));
      gamcen += -0.5 * sum(nobscen(arma::find(nobscen==i)))/i * val;
    }
  }

  double logpost = gamobs + gamcen + cond + obs + cen - prior;
  List Result = List::create(Named("Gamma")=Gamma, _["logpost"]=logpost);

  return Result;
}

//[[Rcpp::export]]
double ar1postdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, arma::vec & r12, double & ar11, double & ar12, double & ar1, arma::vec & alpha, List & idx, double & ar1priorBETA, arma::uvec & cenIM, arma::uvec & obsIM, double & Gammastr_ar1) {

  int m = M*J;
  int mm = M+1 ; int jj = J-1;

  arma::mat Sigma; 
  if (Gammastr_ar1==3) {
    Sigma = AR1cal(ar1,m,Gammastr_ar1);
  }
  else if (Gammastr_ar1==4) {
    double ar1double = r12(0);
    Sigma = kron(AR1cal(ar1double,M,Gammastr_ar1),AR1cal(ar1,J,Gammastr_ar1)); 
  } else {
    Sigma = kron(AR1cal(ar1,M,Gammastr_ar1),Gammacal(r12,jj)); 
  }
  arma::mat Sigmainv = inv(Sigma);
  arma::mat Gammainv(m+1,m+1);
  Gammainv(0,0) = as_scalar(alpha.t() * Sigmainv * alpha) + 1;
  Gammainv(arma::span(1,m),arma::span(1,m)) = Sigmainv;
  Gammainv(0,arma::span(1,m)) = - alpha.t() * Sigmainv ; Gammainv(arma::span(1,m),0) = - Sigmainv * alpha;
  
  arma::mat Gamma = inv(Gammainv);

  arma::vec nobsobs = idx["nobsobs"] ; arma::vec nobscen = idx["nobscen"]; 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double prior = (ar11-1)*log(ar1) + (ar12-1)*log(1-ar1);

  double gamobs = 0.0; double gamcen = 0.0;
  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<obsIM.n_elem; p++) { 
        arma::uvec ttmp = find_finite(Xobs(obsIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamobs += -0.5 * val;
      }
    } else {
      double val; double sign;  
      log_det(val, sign, Gamma(arma::span(0,i*J),arma::span(0,i*J)));
      gamobs += -0.5 * sum(nobsobs(arma::find(nobsobs==i)))/i * val;  
    }
  }

  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<cenIM.n_elem; p++) { 
        arma::uvec ttmp = 1+find_finite(Xcenm(cenIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamcen += -0.5 * val;
      }
    } else {
      double val; double sign;
      log_det(val, sign, Gamma(arma::span(1,i*J),arma::span(1,i*J)));
      gamcen += -0.5 * sum(nobscen(arma::find(nobscen==i)))/i * val;
    }
  }

  double logpost;
  if (ar1priorBETA==1.0)
    logpost = gamobs + gamcen + cond + obs + cen + prior;
  else 
    logpost = gamobs + gamcen + cond + obs + cen;

  return logpost;
}

//[[Rcpp::export]]
double alphapostdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & M, int & J, arma::mat & Xobs, arma::mat & Xcen, arma::mat & Xcenm, arma::vec & Xcent, arma::vec & r12, double & ar1, double & alphaa, arma::vec & alpha_old, int & ind, arma::vec & alpha0, arma::mat & Sigalpha0, double & mh, List & idx, arma::uvec & cenIM, arma::uvec & obsIM , double & Gammastr_ar1) {

  int m = M*J;
  arma::vec alpha = alpha_old;
  if (mh==0.0)
    alpha(ind-1) = alphaa;

  arma::mat Sigma; int jj = J-1;
  if (Gammastr_ar1==3) {
    Sigma = AR1cal(ar1,m,Gammastr_ar1);
  }
  else if (Gammastr_ar1==4) {
    double r12double = r12(0);
    Sigma = kron(AR1cal(ar1,M,Gammastr_ar1),AR1cal(r12double,J,Gammastr_ar1)); 
  } else {
    Sigma = kron(AR1cal(ar1,M,Gammastr_ar1),Gammacal(r12,jj)); 
  }
  arma::mat Sigmainv = inv(Sigma);
  arma::mat Gammainv(m+1,m+1);
  Gammainv(0,0) = as_scalar(alpha.t() * Sigmainv * alpha) + 1;
  Gammainv(arma::span(1,m),arma::span(1,m)) = Sigmainv;
  Gammainv(0,arma::span(1,m)) = - alpha.t() * Sigmainv ; Gammainv(arma::span(1,m),0) = - Sigmainv * alpha;
  
  arma::mat Gamma = inv(Gammainv);

  arma::vec nobsobs = idx["nobsobs"] ; arma::vec nobscen = idx["nobscen"]; 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double prior;
  if (mh==0.0)
    prior = 0.5 * (1/Sigalpha0(ind-1,ind-1)) * (alphaa-alpha0(ind-1))*(alphaa-alpha0(ind-1));
  else 
    prior = 0.5 * as_scalar((alpha-alpha0).t() * inv(Sigalpha0) * (alpha-alpha0));

  double gamobs = 0.0; double gamcen = 0.0;
  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<obsIM.n_elem; p++) { 
        arma::uvec ttmp = find_finite(Xobs(obsIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamobs += -0.5 * val;
      }
    } else {
      double val; double sign;  
      log_det(val, sign, Gamma(arma::span(0,i*J),arma::span(0,i*J)));
      gamobs += -0.5 * sum(nobsobs(arma::find(nobsobs==i)))/i * val;  
    }
  }

  for (int i=1; i<M+1; i++)
  {
    if (i==M) {
      for (int p=0; p<cenIM.n_elem; p++) { 
        arma::uvec ttmp = 1+find_finite(Xcenm(cenIM(p),arma::span::all));

        double val; double sign;
        log_det(val, sign, Gamma.submat(ttmp,ttmp));
        gamcen += -0.5 * val;
      }
    } else {
      double val; double sign;
      log_det(val, sign, Gamma(arma::span(1,i*J),arma::span(1,i*J)));
      gamcen += -0.5 * sum(nobscen(arma::find(nobscen==i)))/i * val;
    }
  }

  double logpost = gamobs + gamcen + cond + obs + cen - prior;
  return logpost;
}

//[[Rcpp::export]]
double mupostdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & muu, arma::vec & mu_old, int & ind, arma::vec & mu0, arma::mat & Sig0, arma::vec & sig2, arma::mat & Gamma, arma::cube & W, double & mh, List & idx) {
  arma::vec mu = mu_old;
  if (mh==0.0)
  	mu(ind-1) = muu;
  arma::mat ymu(n,m);
  if (W.size()>1) {
    arma::mat Wmu(m,n);
    for (int i=0; i<n; i++)
      Wmu(arma::span::all,i) = W.slice(i) * mu;
    ymu = Y - Wmu.t();
  } else {
    arma::mat Yt = Y.t();
    arma::mat ymut(m,n);
    for (int i=0; i<n; i++)
      ymut(arma::span::all,i) = Yt(arma::span::all,i) - mu;
    ymu = ymut.t();
  }
  arma::mat SSIG(n,m); 
  if (sig2.size()==1) {
    double ss = sig2(0);
    SSIG.fill(ss);
  } else {
    arma::mat Sig2(n,J); Sig2.each_row() = reshape(sig2,1,J);
    for (int i=0;i<M;i++)
      SSIG(arma::span::all,arma::span(J*i,J*(i+1)-1)) = Sig2;  
  }
  arma::mat Xm = ymu / sqrt(SSIG);
  for (int i=0; i<m; i++)
    Xm(arma::span::all,i) = Xm(arma::span::all,i) * sqrt(Gamma(i+1,i+1));

  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::mat Xcen  = X.rows(arma::find(v==0));
  arma::mat Xcenm = Xm.rows(arma::find(v==0));
  arma::vec Xcent = Xt.elem(arma::find(v==0));

  arma::mat ymu2 = Xm % Xm;
  double datalik  = -0.5 * accu(ymu2.elem(find_finite(ymu2)));
  double prior = -0.5 * accu((mu-mu0).t() * inv(Sig0) * (mu-mu0));
  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);
  double logpost = obs + cen + cond + datalik + prior;
  return logpost;
}

//[[Rcpp::export]]
double sig2postdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & sig2, double & c0, double & d0, arma::vec & mu, arma::mat & Gamma, arma::cube & W, List & idx) {

  arma::mat ymu(n,m);
  if (W.size()>1) {
    arma::mat Wmu(m,n);
    for (int i=0; i<n; i++)
      Wmu(arma::span::all,i) = W.slice(i) * mu;
    ymu = Y - Wmu.t();
  } else {
    arma::mat Yt = Y.t();
    arma::mat ymut(m,n);
    for (int i=0; i<n; i++)
      ymut(arma::span::all,i) = Yt(arma::span::all,i) - mu;
    ymu = ymut.t();
  }

  arma::mat Xm = ymu / sqrt(sig2);
  for (int i=0; i<m; i++)
    Xm(arma::span::all,i) = Xm(arma::span::all,i) * sqrt(Gamma(i+1,i+1));

  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::mat Xcen  = X.rows(arma::find(v==0));
  arma::mat Xcenm = Xm.rows(arma::find(v==0));
  arma::vec Xcent = Xt.elem(arma::find(v==0));

  arma::mat ymu2 = ymu % ymu;

  int ii = 0;
  for (int j=0 ; j<n ; j++)
  {
  	for (int i = 0 ; i < m ; i++) 
  	{
  		ii += NumericVector::is_na(Y(j,i));
  	}
  }

  double shape = c0 + 0.5 * (n*m - ii);
  double rate = 0 ; rate  = d0 + 0.5 * accu(ymu2.elem(find_finite(ymu2))); 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double logpost = - rate/sig2 - (shape + 1)*log(sig2) + cond + obs + cen;

  return logpost;
}

//[[Rcpp::export]]
double sig2multpostdist_cpp(arma::mat & Y, arma::vec & v, int & n, int & m, int & M, int & J, arma::vec & Xt, double & sig22, arma::vec & sig2_old, int & ind, double & c0, double & d0, arma::vec & mu, arma::mat & Gamma, arma::cube & W, List & idx) {

  arma::vec sig2 = sig2_old;
  sig2(ind-1) = sig22;

  arma::mat ymu(n,m);
  if (W.size()>1) {
    arma::mat Wmu(m,n);
    for (int i=0; i<n; i++)
      Wmu(arma::span::all,i) = W.slice(i) * mu;
    ymu = Y - Wmu.t();
  } else {
    arma::mat Yt = Y.t();
    arma::mat ymut(m,n);
    for (int i=0; i<n; i++)
      ymut(arma::span::all,i) = Yt(arma::span::all,i) - mu;
    ymu = ymut.t();
  }

  arma::mat SSIG(n,m); 
  if (sig2.size()==1) {
    double ss = sig2(0);
    SSIG.fill(ss);
  } else {
    arma::mat Sig2(n,J); Sig2.each_row() = reshape(sig2,1,J);
    for (int i=0;i<M;i++)
      SSIG(arma::span::all,arma::span(J*i,J*(i+1)-1)) = Sig2;  
  }

  arma::mat Xm = ymu / sqrt(SSIG);
  for (int i=0; i<m; i++)
    Xm(arma::span::all,i) = Xm(arma::span::all,i) * sqrt(Gamma(i+1,i+1));

  
  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::mat Xcen  = X.rows(arma::find(v==0));
  arma::mat Xcenm = Xm.rows(arma::find(v==0));
  arma::vec Xcent = Xt.elem(arma::find(v==0));

  arma::mat ymu2 = ymu % ymu;

  int ii = 0;
  for (int j=0 ; j<n ; j++)
  {
    for (int i = 0 ; i < m ; i++) 
    {
      ii += NumericVector::is_na(Y(j,i));
    }
  }

  double shape = c0 + 0.5 * (n*m - ii);
  double rate = 0 ; rate  = d0 + 0.5 * accu(ymu2.elem(find_finite(ymu2))); 

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);
  double cen  = censliklog_cpp(Xcenm,m,M,J,Gamma, idx);

  double logpost = - rate/sig22 - (shape + 1)*log(sig22) + cond + obs + cen;

  return logpost;
}

//[[Rcpp::export]]
double lambdapostdist_cpp(arma::vec & Time, arma::mat & Z, arma::vec & v, arma::vec & delta, int & n, int & m, int & M, int & J, double & Xtmax, double & Xtmin, arma::mat & Xm, arma::mat & Xcenm, double & lambdaa, arma::vec & lambda_old, int & ind, arma::vec & timeint, double & a0, double & b0, arma::vec & beta, arma::mat & Gamma, List & idx) {

  arma::vec cumsum = arma::zeros(n); arma::vec survind = arma::zeros(n); arma::vec surv = arma::zeros(n); 
  double failtime=0.0 ; double failtimeobs=0.0; double shape=0.0; double rate=0.0; double lambrelated=0.0; 
  double cond,obs,logpost ; cond=obs=logpost=0.0; 

  arma::vec lambda = lambda_old;
  
  if (timeint.size()>1) {
    lambda(ind-1) = lambdaa;

    for (int i=0; i<n; i++)
    {
      if (delta(i)==1) {
        cumsum(i) = 0;
      } else {
        cumsum(i) = sum(lambda(arma::span(0,delta(i)-2)) % (timeint(arma::span(1,delta(i)-1)) - timeint(arma::span(0,delta(i)-2))));
      }
      survind(i) = -( lambda(delta(i)-1) * (Time(i)-timeint(delta(i)-1)) + cumsum(i) ) * exp(accu(Z(i,arma::span::all) % beta.t()));
    }
    surv = exp(survind);
    failtime = accu(survind);
    failtimeobs = accu(survind.elem(arma::find(v==1)));

    shape = a0 + accu(v.elem(arma::find(delta==ind)));
    rate = b0;
    lambrelated = (shape-1) * log(lambdaa) - lambdaa * rate;

  } else {
    arma::vec Tobs = Time(arma::find(v==1));
    arma::mat Zobs = Z.rows(arma::find(v==1));
    arma::vec ezb(Zobs.n_rows);
    
    shape = a0 + accu(v);
    ezb = exp(Zobs * beta);
    rate = b0 + accu(Tobs % ezb);
    failtime = 0.0;
    lambrelated = (shape-1) * log(lambdaa) - lambdaa * rate;
    surv = exp(-lambda * Time % exp(Z * beta));
  }

  arma::vec Xt = arma::zeros(n); 
  for (int i=0; i<n; i++) {
  	Xt(i) = R::qnorm(1-surv(i), 0.0, 1.0, TRUE, FALSE);
  	if (!arma::is_finite(Xt(i))){
      if (surv(i)<0.1)
        Xt(i) = Xtmax;
      else 
        Xt(i) = Xtmin;
      //Xt(i) = 8.2;
      //Rcout << "Xt==8.2\n";
      //Rcout << Xt(i) << "\n";
    }
  }

  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::vec Xcent = Xt.elem(arma::find(v==0));

  cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);

  logpost = lambrelated + cond + obs + failtimeobs;

  return logpost;
}

//[[Rcpp::export]]
double betapostdist_cpp(arma::vec & Time, arma::mat & Z, arma::vec & v, arma::vec & delta, int & n, int & m, int & M, int & J, double & Xtmax, double & Xtmin,  arma::mat & Xm, arma::mat & Xcenm, double & betaa, arma::vec & beta_old, int & ind, arma::vec & timeint, arma::vec & beta0, arma::mat & Sigbeta0, arma::vec & lambda, arma::mat & Gamma, List & idx) {

  arma::vec beta = beta_old; beta(ind-1) = betaa;
  
  arma::vec cumsum(n);
  arma::vec survind(n);
  arma::vec surv(n); 
  arma::mat Zobs = Z.rows(arma::find(v==1));
  double failtime=0.0; 
  
  if (timeint.size()>1) {

    for (int i=0; i<n; i++)
    {
      if (delta(i)==1) {
        cumsum(i) = 0;
      } else {
        cumsum(i) = sum(lambda(arma::span(0,delta(i)-2)) % (timeint(arma::span(1,delta(i)-1)) - timeint(arma::span(0,delta(i)-2))));
      }
      survind(i) = -( lambda(delta(i)-1) * (Time(i)-timeint(delta(i)-1)) + cumsum(i) ) * exp(accu(Z(i,arma::span::all) % beta.t()));
    }
    surv = exp(survind);
    failtime = accu(Zobs * beta) + accu(survind.elem(arma::find(v==1)));
  
  } else {
    arma::vec Tobs = Time(arma::find(v==1));
    arma::vec ezb(Zobs.n_rows); arma::vec zb(Zobs.n_rows);
    
    zb = Zobs * beta;
    ezb = exp(zb);
    failtime = 0.0; // need to change so that lambda is a vector
    surv = exp(-lambda * Time % exp(Z * beta));
  }

  double prior = -0.5 * accu((beta-beta0).t() * inv(Sigbeta0) * (beta-beta0));

  arma::vec Xt(n);
  for (int i=0; i<n; i++) {
    double xx = R::qnorm(1-surv(i), 0.0, 1.0, TRUE, FALSE);
    if (!arma::is_finite(xx)) {
      if (surv(i)<0.1)
        Xt(i) = Xtmax;
      else 
        Xt(i) = Xtmin;      
    } else {
      Xt(i) = xx;
    } 
      
  }

  arma::mat X(n,m+1);
  X(arma::span::all,0) = Xt;
  X(arma::span::all,arma::span(1,m)) = Xm;

  arma::mat Xobs  = X.rows(arma::find(v==1));
  arma::vec Xcent = Xt.elem(arma::find(v==0));

  double cond = survcdflog_cpp(Xcenm,Xcent,m,n,M,J,Gamma,idx);
  double obs  = obsliklog_cpp(Xobs, m, M,J, Gamma, idx);

  double logpost = failtime + prior + cond + obs;

  return logpost;
}


//[[Rcpp::export]]
arma::vec slice_sample_cpp(String & param_type, arma::vec & x0, double & w, double & slice_m, double & lower, double & upper, List & params) { 
  
  double J, K;
  arma::vec L(x0.size()), R(x0.size()), param_new(x0.size());
  arma::vec paramvec = param_new = x0;

  List idx = params["idx"]; 

  for (int j=0; j < x0.size(); j++) {

    double logy=0.0; double logz=0.0;

    double param_old = paramvec(j); 
    int jj = j + 1;

    // Calculate initial horizontal interval;
    L(j) = param_old - R::runif(0.0, w);
    R(j) = L(j) + w;
    // Truncate bounds to support of the parameter space;
    L(j) = std::max(L(j),lower);
    R(j) = std::min(R(j),upper);
    // Step out;
    J = floor(slice_m * R::runif(0.0,1.0));
    K = (slice_m-1)-J;

    double funL, funR;
    if (param_type=="mu") {
      arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
      arma::vec mu0 = params["mu0"] ; arma::mat Sig0 = params["Sig0"] ; arma::vec sig2 = params["sig2"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ; double mh=0.0;
    
      logy = mupostdist_cpp(Y, v,n,m,M,J, Xt, param_old, paramvec, jj, mu0,Sig0,sig2,Gamma,W,mh, idx);
      funL = mupostdist_cpp(Y, v,n,m,M,J, Xt, L(j), paramvec, jj, mu0,Sig0,sig2,Gamma,W,mh, idx);
      funR = mupostdist_cpp(Y, v,n,m,M,J, Xt, R(j), paramvec, jj, mu0,Sig0,sig2,Gamma,W,mh, idx);
    } else if (param_type=="lambda") {
      arma::vec Time = params["Time"] ; arma::mat Z=params["Z"] ; arma::vec v=params["v"] ; arma::vec delta=params["delta"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; 
      double Xtmax=params["Xtmax"] ; double Xtmin=params["Xtmin"];
      arma::mat Xm=params["Xm"] ; arma::mat Xcenm=params["Xcenm"] ; arma::mat Gamma = params["Gamma"] ; 
      arma::vec timeint=params["timeint"] ; double a0=params["a0"] ; double b0=params["b0"] ; arma::vec beta=params["beta"] ; 

      logy = lambdapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, param_old, paramvec, jj, timeint,a0,b0,beta,Gamma, idx);
      funL = lambdapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, L(j), paramvec, jj, timeint,a0,b0,beta,Gamma, idx);
      funR = lambdapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, R(j), paramvec, jj, timeint,a0,b0,beta,Gamma, idx);
    } else if (param_type=="Gamma") {
      arma::mat Y = params["Y"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
      double sigr2=params["sigr2"] ; double mh=0.0;
      arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"];

      logy = r12postdist_cpp(Y,n,m,M,J, Xobs,Xcen,Xcenm,Xcent, param_old, paramvec, jj, sigr2, mh, idx,cenIM,obsIM)["logpost"];
      funL = r12postdist_cpp(Y,n,m,M,J, Xobs,Xcen,Xcenm,Xcent, L(j), paramvec, jj, sigr2, mh, idx,cenIM,obsIM)["logpost"];
      funR = r12postdist_cpp(Y,n,m,M,J, Xobs,Xcen,Xcenm,Xcent, R(j), paramvec, jj, sigr2, mh, idx,cenIM,obsIM)["logpost"];
    } else if (param_type=="Gammastr") {
      arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
      arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
      double sigr2=params["sigr2"] ; double ar1=params["ar1"] ; arma::vec alpha=params["alpha"] ; double mh=0.0;
      arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"]; double Gammastr_ar1=params["Gammastr_ar1"];

      logy = Gammapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, param_old, paramvec, jj, sigr2, ar1, alpha, mh, idx,cenIM,obsIM,Gammastr_ar1)["logpost"];
      funL = Gammapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, L(j), paramvec, jj, sigr2, ar1, alpha, mh, idx,cenIM,obsIM,Gammastr_ar1)["logpost"];
      funR = Gammapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, R(j), paramvec, jj, sigr2, ar1, alpha, mh, idx,cenIM,obsIM,Gammastr_ar1)["logpost"];
    } else if (param_type=="ar1") {
      arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
      arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
      arma::vec r12=params["r12"] ; arma::vec alpha=params["alpha"] ; 
      double ar11=params["ar11"] ; double ar12=params["ar12"] ; 
      arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"]; double Gammastr_ar1=params["Gammastr_ar1"];
      double ar1priorBETA = params["ar1priorBETA"] ;

      logy = ar1postdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar11, ar12, param_old, alpha, idx, ar1priorBETA,cenIM,obsIM,Gammastr_ar1);
      funL = ar1postdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar11, ar12, L(j), alpha, idx, ar1priorBETA,cenIM,obsIM,Gammastr_ar1);
      funR = ar1postdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar11, ar12, R(j), alpha, idx, ar1priorBETA,cenIM,obsIM,Gammastr_ar1);
    } else if (param_type=="alpha") {
      arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
      arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
      arma::vec r12=params["r12"] ; double ar1=params["ar1"] ; arma::vec alpha0=params["alpha0"] ; arma::mat Sigalpha0=params["Sigalpha0"] ; double mh=0.0;
      arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"]; double Gammastr_ar1=params["Gammastr_ar1"];

      logy = alphapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar1, param_old, paramvec, jj, alpha0, Sigalpha0, mh, idx,cenIM,obsIM,Gammastr_ar1);
      funL = alphapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar1, L(j), paramvec, jj, alpha0, Sigalpha0, mh, idx,cenIM,obsIM,Gammastr_ar1);
      funR = alphapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar1, R(j), paramvec, jj, alpha0, Sigalpha0, mh, idx,cenIM,obsIM,Gammastr_ar1);
    } else if (param_type=="sig2") {
      arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
      double c0=params["c0"] ; double d0=params["d0"] ; arma::vec mu=params["mu"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ;
    
      logy = sig2postdist_cpp(Y,v,n,m,M,J, Xt, param_old, c0,d0,mu,Gamma,W, idx);
      funL = sig2postdist_cpp(Y,v,n,m,M,J, Xt, L(j), c0,d0,mu,Gamma,W, idx);
      funR = sig2postdist_cpp(Y,v,n,m,M,J, Xt, R(j), c0,d0,mu,Gamma,W, idx);
    } else if (param_type=="sig2mult") {
      arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
      double c0=params["c0"] ; double d0=params["d0"] ; arma::vec mu=params["mu"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ;
    
      logy = sig2multpostdist_cpp(Y,v,n,m,M,J, Xt, param_old, paramvec, jj, c0,d0,mu,Gamma,W, idx);
      funL = sig2multpostdist_cpp(Y,v,n,m,M,J, Xt, L(j), paramvec, jj, c0,d0,mu,Gamma,W, idx);
      funR = sig2multpostdist_cpp(Y,v,n,m,M,J, Xt, R(j), paramvec, jj, c0,d0,mu,Gamma,W, idx);
    } else if (param_type=="beta") {
      arma::vec Time = params["Time"] ; arma::mat Z=params["Z"] ; arma::vec v=params["v"] ; arma::vec delta=params["delta"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; 
      double Xtmax=params["Xtmax"] ; double Xtmin=params["Xtmin"];
      arma::mat Xm=params["Xm"] ; arma::mat Xcenm=params["Xcenm"] ; arma::mat Gamma = params["Gamma"] ; 
      arma::vec timeint=params["timeint"] ; arma::vec beta0=params["beta0"] ; arma::mat Sigbeta0=params["Sigbeta0"] ; arma::vec lambda=params["lambda"] ; 

      logy = betapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, param_old, paramvec, jj, timeint,beta0,Sigbeta0,lambda,Gamma, idx);
      funL = betapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, L(j), paramvec, jj, timeint,beta0,Sigbeta0,lambda,Gamma, idx);
      funR = betapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, R(j), paramvec, jj, timeint,beta0,Sigbeta0,lambda,Gamma, idx);
    } else if (param_type=="missing") {
      arma::mat Yorigin=params["Yorigin"]; arma::vec missind=params["missind"]; arma::vec v=params["v"]; arma::vec Xt=params["Xt"]; arma::vec Xcent=params["Xcent"]; arma::cube W=params["W"]; int n=params["n"]; int m=params["m"]; int M=params["M"]; int J=params["J"];
      arma::mat Gamma=params["Gamma"]; arma::vec mu=params["mu"]; arma::vec sig2=params["sig2"]; 
     
      logy = missingliklog_cpp(Yorigin, param_old, jj, missind, paramvec, v, Xt, Xcent, W, n, m,M,J, Gamma, mu, sig2, idx);
      funL = missingliklog_cpp(Yorigin, L(j), jj, missind, paramvec, v, Xt, Xcent, W, n, m,M,J, Gamma, mu, sig2, idx);
      funR = missingliklog_cpp(Yorigin, R(j), jj, missind, paramvec, v, Xt, Xcent, W, n, m,M,J, Gamma, mu, sig2, idx);
     
     } 

    // draw uniformly from [0, y]
    logz = logy - R::rexp(1.0);
      
    while ( J > 0 && L(j) > lower && funL > logz )
    {
      L(j) = L(j) - w; 
      if (L(j) <= lower) 
        L(j) = lower;
      J = J-1;
    }
    while ( K > 0 && R(j) < upper && funR > logz )
    {
      R(j) = R(j) + w;
      if (R(j) >= upper) 
        R(j) = upper;
      K = K-1;
    }

    // shrinkage procedure
    int cnt = 0;
    do {
      cnt++;
      paramvec(j) = R::runif(L(j),R(j));

      double funnew;
      if (param_type=="mu") {
        arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
        arma::vec mu0 = params["mu0"] ; arma::mat Sig0 = params["Sig0"] ; arma::vec sig2 = params["sig2"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ; double mh=0.0;
      
        funnew = mupostdist_cpp(Y, v,n,m,M,J, Xt, paramvec(j), paramvec, jj, mu0,Sig0,sig2,Gamma,W,mh, idx);
      } else if (param_type=="lambda") {
        arma::vec Time = params["Time"] ; arma::mat Z=params["Z"] ; arma::vec v=params["v"] ; arma::vec delta=params["delta"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; 
        double Xtmax=params["Xtmax"] ; double Xtmin=params["Xtmin"];
        arma::mat Xm=params["Xm"] ; arma::mat Xcenm=params["Xcenm"] ; arma::mat Gamma = params["Gamma"] ; 
        arma::vec timeint=params["timeint"] ; double a0=params["a0"] ; double b0=params["b0"] ; arma::vec beta=params["beta"] ; 

        funnew = lambdapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, paramvec(j), paramvec, jj, timeint,a0,b0,beta,Gamma, idx);
      } else if (param_type=="Gamma") {
        arma::mat Y = params["Y"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
        double sigr2=params["sigr2"] ; double mh=0.0;
        arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"];

        funnew = r12postdist_cpp(Y,n,m,M,J, Xobs,Xcen,Xcenm,Xcent, paramvec(j), paramvec, jj, sigr2, mh, idx,cenIM,obsIM)["logpost"];
      } else if (param_type=="Gammastr") {
        arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
        arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
        double sigr2=params["sigr2"] ; double ar1=params["ar1"] ; arma::vec alpha=params["alpha"] ; double mh=0.0;
        arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"]; double Gammastr_ar1=params["Gammastr_ar1"];

        funnew = Gammapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, paramvec(j), paramvec, jj, sigr2, ar1, alpha, mh, idx,cenIM,obsIM,Gammastr_ar1)["logpost"];
      } else if (param_type=="ar1") {
        arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
        arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
        arma::vec r12=params["r12"] ; arma::vec alpha=params["alpha"] ; 
        double ar11=params["ar11"] ; double ar12=params["ar12"] ; 
        arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"];  double Gammastr_ar1=params["Gammastr_ar1"];
        double ar1priorBETA=params["ar1priorBETA"] ;

        funnew = ar1postdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar11, ar12, paramvec(j), alpha, idx, ar1priorBETA,cenIM,obsIM,Gammastr_ar1);
      } else if (param_type=="alpha") {
        arma::mat Y = params["Y"] ; arma::vec v=params["v"] ; int n = params["n"] ; int M=params["M"]; int J=params["J"]; 
        arma::mat Xobs=params["Xobs"] ; arma::mat Xcen=params["Xcen"] ; arma::mat Xcenm=params["Xcenm"] ; arma::vec Xcent = params["Xcent"] ;
        arma::vec r12=params["r12"] ; double ar1=params["ar1"] ; arma::vec alpha0=params["alpha0"] ; arma::mat Sigalpha0=params["Sigalpha0"] ; double mh=0.0;
        arma::uvec cenIM=params["cenIM"]; arma::uvec obsIM=params["obsIM"]; double Gammastr_ar1=params["Gammastr_ar1"];

        funnew = alphapostdist_cpp(Y,v,n,M,J, Xobs,Xcen,Xcenm,Xcent, r12, ar1, paramvec(j), paramvec, jj, alpha0, Sigalpha0, mh, idx,cenIM,obsIM,Gammastr_ar1);
      } else if (param_type=="sig2") {
        arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
        double c0=params["c0"] ; double d0=params["d0"] ; arma::vec mu=params["mu"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ;
      
        funnew = sig2postdist_cpp(Y,v,n,m,M,J, Xt, paramvec(j), c0,d0,mu,Gamma,W, idx);
      } else if (param_type=="sig2mult") {
        arma::mat Y = params["Y"] ; arma::vec v = params["v"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; arma::vec Xt = params["Xt"] ;
        double c0=params["c0"] ; double d0=params["d0"] ; arma::vec mu=params["mu"] ; arma::mat Gamma = params["Gamma"] ; arma::cube W = params["W"] ;
      
        funnew = sig2multpostdist_cpp(Y,v,n,m,M,J, Xt, paramvec(j), paramvec, jj, c0,d0,mu,Gamma,W, idx);
      } else if (param_type=="beta") {
        arma::vec Time = params["Time"] ; arma::mat Z=params["Z"] ; arma::vec v=params["v"] ; arma::vec delta=params["delta"] ; int n = params["n"] ; int m = params["m"] ; int M=params["M"]; int J=params["J"]; 
        double Xtmax=params["Xtmax"] ; double Xtmin=params["Xtmin"];
        arma::mat Xm=params["Xm"] ; arma::mat Xcenm=params["Xcenm"] ; arma::mat Gamma = params["Gamma"] ; 
        arma::vec timeint=params["timeint"] ; arma::vec beta0=params["beta0"] ; arma::mat Sigbeta0=params["Sigbeta0"] ; arma::vec lambda=params["lambda"] ; 

        funnew = betapostdist_cpp(Time,Z,v,delta,n,m,M,J, Xtmax,Xtmin, Xm,Xcenm, paramvec(j), paramvec, jj, timeint,beta0,Sigbeta0,lambda,Gamma, idx);
      } else if (param_type=="missing") {
        arma::mat Yorigin=params["Yorigin"]; arma::vec missind=params["missind"]; arma::vec v=params["v"]; arma::vec Xt=params["Xt"]; arma::vec Xcent=params["Xcent"]; arma::cube W=params["W"]; int n=params["n"]; int m=params["m"]; int M=params["M"]; int J=params["J"];
        arma::mat Gamma=params["Gamma"]; arma::vec mu=params["mu"]; arma::vec sig2=params["sig2"]; 
       
        funnew = missingliklog_cpp(Yorigin, paramvec(j), jj, missind, paramvec, v, Xt, Xcent, W, n, m, M,J, Gamma, mu, sig2, idx);
      } 
      
      if ( funnew > logz )
        break;
      if ( paramvec(j) < param_old )
        L(j) = paramvec(j);
      else
        R(j) = paramvec(j);
    } while (cnt<1e4);
      if (cnt==1e4) ::Rf_error("slice_sample_cpp loop did not finish");

      if (-0.0000000001 <= L(j) - R(j) && L(j) - R(j) <= 0.0000000001)
        paramvec(j) = 0.5 * (L(j) + R(j));
    }

  return paramvec;
}
