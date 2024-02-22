# rm(list=ls())

# set.seed(1)

## Code used by UNC computing cluster to pull in inputs from Bash shell;
options(echo=TRUE)
args <- commandArgs(trailingOnly = TRUE)
node.idx = as.numeric(args[1]);
node.idx <- Sys.getenv("SLURM_ARRAY_TASK_ID")

## set node.idx for running in windows environment (i.e., not on UNC computing cluster);
if (.Platform$OS.type == "windows") { node.idx = 1 }
# if (!is.na(node.idx)) node.idx = formatC(node.idx, width = 4, format = "d", flag = "0")

node.idx = ifelse(is.na(node.idx), 3, node.idx)
nchain = node.idx %% 5
nnn = ceiling(node.idx/5) ; nnn = nnn%%3

cat('node.idx : ',node.idx,"\n")
cat('nchain : ',nchain,"\n")
cat('nnn : ',nnn,"\n")

library(mvtnorm)
library(lattice)
library(data.table)
library(Rcpp)
library(mice) # MISSING AT RADOM
library(gee) ; library(geepack) ; # library(mmm2)
library(tidyr)

## ----------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------------------------------- ##

##### priors : hyper-parameters ##### 
### a0, b0 : lambda
### mu0, Sig0 : mu
### c0, d0 : sig2
### beta0, Sigbeta0 : beta
### sigr2 : r12
### ar11, ar12, ar1priorBETA : rho (structured)
### alpha0, Sigalpha0 : alpha (structured)
##### ------------------------ #####


## ----------------------------------------------------------------------------------------- ##
## ----------------------------------------------------------------------------------------- ##
structured = TRUE
MLEestimate=TRUE

setwd('where you have your functions')
true = readRDS(paste0("truevalues.RData"))

cpp=TRUE ; MAR=FALSE
Gammastr_ar1=4
if (cpp){
  Sys.setenv("PKG_CXXFLAGS"="-std=c++11")
  sourceCpp('./Cpp/MCMC_cppMAR_CPPWRAP.cpp')
  sourceCpp('./Cpp/Gcopula_jointMAR.cpp')
}

###------------###
###--- DATA ---###
###------------###
n = true$dim$n ; n1 = true$dim$n1; n0 = true$dim$n0 ; 
M = true$dim$m ; # number of longitudinal measures
J = true$dim$J ; # number of time points
K = true$dim$K ; q = true$dim$q ; zm = true$dim$zm ;

true$structured = structured

m = M*J

Posbb = matrix(c(1,0.3,0.3,1),2,2) ; Negbb = matrix(c(1,-0.3,-0.3,1),2,2)
bb = array(NA,dim=c(n,2,M))
for (i in 1:M) {
  if (i%%2==1) {
    bb[,,i] = rmvnorm(n,numeric(2),Posbb)  
  } else {
    bb[,,i] = rmvnorm(n,numeric(2),Negbb)
  }
}
ttij = seq(0.5, 1,length=J) ; tij = rep(ttij,each=M)
WW = matrix(rnorm(n*(q-2)),nr=n)

## Covariate
Z = matrix(0,n,zm)
Z = matrix(rnorm(n*zm),n,zm)
betat = true$true$betat

W = array(NA, dim=c(m,q*M,n))
timeidx = (1:(q*M))[(1:(q*M) %% q ==2)]

invisible(sapply(1:n, function(i) {
  tmp = matrix(0,m,q*M)
  for (j in 1:M) {
    tmp[seq(j,m,by=M),timeidx[j]] = ttij # time covariate
    if (q>2){ # standard normal covariates
      for (k in 1:(q-2)) {
        tmp[seq(j,m,by=M),timeidx[j]+k] = WW[i,k]
      }
    }
    tmp[seq(j,m,by=M),timeidx[j]-1] = 1 + bb[i,1,j] + bb[i,2,j] * ttij   # intercept  
  }
  W[,,i] <<-  tmp
}))

## Copula
Gammat = true$true$Gammat ; Gammat = cov2cor(Gammat) ; 
ar1 = true$true$ar1 ; r12t = true$true$rhot
alphat = true$true$alphat ; Sigmat = true$true$Sigmat

X = rmvnorm(n,numeric(m+1),Gammat)

## T
lambdat = true$true$lambdat
timeint = true$timeint

## T
pwe=FALSE
T = - log(1-pnorm(X[,1])) / (lambdat[1] * exp(Z %*% betat))
C = runif(n,0,quantile(T,1))

if (pwe){
  Tgenerate = function(u, timeint, lambda, Z, beta) {
    sapply(2:4, function(j) sapply( 1:length(u), function(i) u[i] <= Tcdf(i, timeint[j], timeint, lambda, Z, beta) ))
  }
  
  tmp = Tgenerate(pnorm(X[,1]), timeint, lambdat, Z, betat)
  delta = sapply(1:n, function(i) which(tmp[i,]==1)[1])
  
  tmp = sapply(1:n, function(i) ifelse(delta[i]==1, 0, sum(lambdat[1:(delta[i]-1)] * (timeint[2:delta[i]]-timeint[1:(delta[i]-1)]))))
  T = ( -(log(1-pnorm(X[,1])))/(exp(Z%*%betat)) + lambdat[delta]*timeint[delta] - tmp ) / (lambdat[delta])
  C = runif(n,0,max(T))
  
  if (max(T) > timeint[K]) {
    timeint[K+1] = max(T) + 1  
  } else {
    timeint = timeint[-(K+1)]
    K = K - 1
    lambdat = lambdat[1:K]
  }
}


## Y 
mut = true$true$mut ; sig2t = true$true$sig2t
etat = true$true$etat

Y = sqrt(sig2t) * X[,2:(m+1)] 
Y = sweep(Y,MARGIN=2, sqrt(diag(Gammat)[-1]),"/")


if (exists('W')) {
  Y = Y + t(apply(W,3, function(x) x %*% etat))
} else {
  Y = sweep(Y, MARGIN = 2, mut,'+')  
}

### MISSING AT RANDOM ###
if (MAR) Y = ampute(Y, mech="MAR")$amp

## Data
dat = data.frame(X, Y=Y, T, C, Tst = pmin(T,C), v=1*(T<C), Z=Z)

timeint = true$timeint = c(0,round(quantile(dat[dat$v==1,'Tst'],c(1/3,2/3)),1),Inf)

dat = data.frame(dat, delta=findInterval(dat$Tst,timeint, left.open=TRUE))
names(dat)[grep('Y',names(dat))] = paste0('Y',rep(1:M,J),rep(1:J,each=M))
names(dat)[grep('X',names(dat))] = c("Xt",paste0('X',rep(1:M,J),rep(1:J,each=M)))

Y = dat[,grep("Y",colnames(dat))]

tmp = cut(dat$Tst,c(seq(from=0, by=1.5, length=J),max(dat$Tst)+10))
ind = rowSums(sapply(1:J, function(i) ifelse(tmp==levels(tmp)[i], i,0)))
#-- INTERMITTENT MISSING : UNBALANCED LONGITUDINAL
tmpind = sample(which(ind==2),size=floor(sum(ind==2)/3)) ; Y[tmpind,c((1+1*M):(2*M))] = NA ; ind[tmpind] = J ; rm(tmpind)
#--
dat$obsind = ind
invisible(sapply(1:(J-1), function(i) Y[ind==i,c((i*M+1):ncol(Y))] <<- NA))

dat[,grep("Y", colnames(dat))] = Y
dat[,grep("X",colnames(dat))][,-1][is.na(Y)] = NA

#-- INTERMITTENT MISSING : UNBALANCED LONGITUDINAL
obsind_idx = which(dat$obsind==J)
tmp = apply(dat[obsind_idx,grep("Y", colnames(dat))],1,is.na) ; ttmp = apply(tmp,2,which)
unique(ttmp) ; 
IMidx = vector('list',length=length(unique(ttmp))) ; names(IMidx) = unique(ttmp)

for (i in 1:length(obsind_idx)) {
  for (j in 1:length(unique(ttmp))) {
    if (identical(unlist(unique(ttmp)[j]),ttmp[[i]])) {
      IMidx[[j]] = c(IMidx[[j]],obsind_idx[i])    
    }
  }
}
IMidx = IMidx[!unlist(lapply(IMidx,is.null))]

IMidx
if (length(IMidx)>1) stop('Intermittent Missing cases more than one')

rm(tmp,ttmp)
tmp = IMidx ; rm(IMidx)
IMidx = vector('list',length=2) ; names(IMidx) = c('colnum','idx')
IMidx$colnum = IMidx$idx = vector('list',length=length(tmp))
aa = names(tmp)

IMidx$colnum = c(as.numeric(substr(aa[aa %like% ":"],1,regexpr(":",aa[aa %like% ":"])[1]-1)):as.numeric(substr(aa[aa %like% ":"],regexpr(":",aa[aa %like% ":"])[1]+1,nchar(aa[aa%like%":"]))))
IMidx$colnum = which(1:(M*J)!=IMidx$colnum)
IMidx$idx = tmp[[1]]
IMidx$fullyobservedidx = obsind_idx[!(obsind_idx %in% IMidx$idx)]
cenidx = which(dat[,'v']==0) ; obsidx = which(dat[,'v']==1)
IMidx$cenidx = which(cenidx %in% IMidx$idx) ; IMidx$obsidx = which(obsidx %in% IMidx$idx)
IMidx$cenfullyobservedidx = which(cenidx %in% IMidx$fullyobservedidx) ; IMidx$obsfullyobservedidx = which(obsidx %in% IMidx$fullyobservedidx)
#--

dat = as.matrix(dat)

for (i in 1:n) {
  W[-(1:(ind[i]*M)),,i] = NA
}

cen = 1 - sum(dat[,'v']==1)/n
cat('censoring proportion :', cen, '\n')

###--------------###
###--- PRIORS ---###
###--------------###
## lambda
# a0 = 1/400  ; b0    = 1/200
a0 = 5 ; b0 = 10
a0 = 10 ; b0 = 20
## mu
if (exists('W')) {
  mu0 = numeric(dim(W)[2])
  Sig0 = 100 * diag(dim(W)[2])
} else {
  W = NA
  mu0 = numeric(m) ; Sig0 = 100 * diag(m)  
}
## sig2
c0 = 3  ; d0    = 20
## Gamma (correlation matrix)
sigr2 = 100
## beta
beta0 = numeric(zm) ; Sigbeta0 = diag(100,zm)
## ar1
ar1priorBETA = 1*FALSE; filename_ar1prior = ifelse(ar1priorBETA,"ar1prior_BETA","ar1prior_Unif")
ar11 = 1/3 ; ar12 = 1/3
## alpha
alpha0 = numeric(m) ; Sigalpha0 = 1 * diag(m)

priors = list(a0=a0, b0=b0, mu0=mu0, Sig0=Sig0, c0=c0, d0=d0, sigr2=sigr2, beta0=beta0, Sigbeta0=Sigbeta0, ar11=ar11, ar12=ar12, alpha0=alpha0, Sigalpha0=Sigalpha0, ar1priorBETA=ar1priorBETA)


saveRDS(list(true=true$true,dat=dat, W=W, priors=priors, timeint=timeint),paste0("./",node.idx,"-Simuldatacpp.RData"))


###------------###
###--- DATA ---###
###------------###
dd = readRDS(paste0("./",node.idx,"-Simuldatacpp.RData"))
dat = dd$dat
W = dd$W

Y = dat[,grep('Y',colnames(dat))]  
m = ncol(Y) ; n = nrow(dat)
X = dat[,grep('X',colnames(dat))] ; Xm = X[,2:(m+1)] ; Xt = X[,1]
corrlength = choose(m,2)
timeint = dd$timeint
K = length(timeint)-1
Z = as.matrix(dat[,grep('Z',colnames(dat))])

###---------------------###
###--- MLE Estimates ---###
###---------------------###

simuldat=readRDS("truevalues202312_FALSE.RData")

## Longitudinal
MLE.estimates = matrix(NA,(q*M+1+2*choose(J,2)) + zm+K, 4) ; colnames(MLE.estimates) = c('Estimate','RobustSE','ll','ul') ; rownames(MLE.estimates) = c(paste0("W",rep(1:M,each=q),rep(1:q,M)),'rho',paste0('r1.',1:choose(J,2)),paste0('r2.',1:choose(J,2)), paste0("Z",1:zm),paste0("lambda",1:K))

if (MLEestimate) {
  
for (i in 1:M) {
  geedat = data.frame(ID = 1:n, Y[,seq(i,length=J,by=M)])
  ww = W[seq(i,length=J,by=M),((i-1)*q+1):(i*q),]
  
  geedatl = gather(geedat,time,y,colnames(geedat)[colnames(geedat) %like% "Y"]) ; geedatl = geedatl[order(geedatl$ID),]
  www = apply(ww,2,c)
  geedatl = data.frame(geedatl, W=www)

  form = as.formula(paste("y ~ ", BBmisc::collapse(paste0("W.",1:q),"+"),"-1"))

  fitun = gee(form, id=ID, data=geedatl, corstr='unstructured') 
  fitstr = geeglm(form, id=ID, data=geedatl, family=gaussian, corstr="ar1")
  fitunglm = geeglm(form, id=ID, data=geedatl, family=gaussian, corstr="unstructured")
  
  MLE.estimates[(1+(i-1)*q):(i*q),1:2] = summary(fitun)$coefficients[,c(1,4)]
  MLE.estimates[(1+(q*M) + 1 + (i-1)*choose(J,2)):(1+(q*M) + (i*choose(J,2))),1] = fitunglm$geese$alpha # fitun$working.correlation
  MLE.estimates[(1+(q*M) + 1 + (i-1)*choose(J,2)):(1+(q*M) + (i*choose(J,2))),2] = sqrt(diag(fitunglm$geese$valpha))
  MLE.estimates['rho',1] = fitstr$geese$alpha
  MLE.estimates['rho',2] = sqrt(fitstr$geese$valpha)
}

## Survival
library(survival) ; library(pch) ; library(eha)

survdat = data.frame(Tst=dat[,'Tst'],v=dat[,'v'],Z=dat[,colnames(dat) %like%'Z'])
form = as.formula(paste("Surv(Tst,event=v) ~ ",BBmisc::collapse(paste0("Z.Z.",1:zm),"+")))
fit = eha::pchreg(form, cuts=timeint, data=survdat)
MLE.estimates[rownames(MLE.estimates) %like% "Z",1] = fit$coefficients 
MLE.estimates[rownames(MLE.estimates) %like% "Z",2] = sqrt(diag(fit$var))
MLE.estimates[rownames(MLE.estimates) %like% "lambda",1] = fit$hazards

}

###----------------###
###--- INITIALS ---###
###----------------###
if (nchain %in% c(0,1)){
  
  sig2   = 1 
  if (!all(is.na(W))) {
    mu = numeric(dim(W)[2])
  } else {
    mu = numeric(m)
  }
  lambda = rep(1,K) 
  beta = numeric(zm) 
  ar1 = 0.5 ; alpha = rep(0.01, m)
  if (true$structured) {
    r12 = numeric(choose(M,2))
    Sigma = ar1_cor(J,ar1) %x% Corrmat(M,r12)$Gamma
    Gammainv = matrix(NA, m+1,m+1) ; Gammainv[1,1] = crossprod(alpha, solve(Sigma,alpha)) + 1 ; Gammainv[-1,-1] = solve(Sigma) ; Gammainv[-1,1] = Gammainv[1,-1] = -solve(Sigma,alpha)
    Gamma = solve(Gammainv)
  } else {
    r12    = numeric(choose(m+1,2)) 
    Gamma = Corrmat(m+1,r12)$Gamma
  }
  
} else if (nchain %in% c(2:4)) {

  sig2   = sd(Y[,1], na.rm=TRUE) 
  if (!all(is.na(W))) {
    mu = colMeans(apply(W,2,colMeans),na.rm=TRUE)
  } else {
    mu = colMeans(Y,na.rm=TRUE)
  }
  lambda = MLE.estimates[rownames(MLE.estimates) %like% "lambda",1] 
  beta = MLE.estimates[rownames(MLE.estimates) %like% "Z",1] 
  ar1 = ifelse(true$structured,abs(cor(X[!(is.na(X[,M+2]) | is.na(X[,2])),2],X[!(is.na(X[,M+2]) | is.na(X[,2])),M+2])),0.5)
  Gamma = cor(X,use="complete.obs")
  Gammainv = solve(Gamma) ; alpha = -solve(Gammainv[2:ncol(Gamma),2:ncol(Gamma)],Gammainv[-1,1])

  if (true$structured) {
    r12 = var(X[2:(M+1),2:(M+1)],na.rm=TRUE)[upper.tri(var(X[2:(M+1),2:(M+1)]))]
    if (is.na(r12)) r12 = rnorm(choose(M,2))
    Sigma = ar1_cor(J,ar1) %x% Corrmat(M,r12)$Gamma
  } else {
    r12    = var(X,na.rm=TRUE)[upper.tri(var(X,na.rm=TRUE))]
    if (is.na(r12)) r12 = rnorm(choose(M,2))
    Gamma = Corrmat(m+1,r12)$Gamma
  }
  
} else {

  sig2   = 1/rgamma(1,c0,rate=d0) 
  if (!all(is.na(W))) {
    mu = rnorm(dim(W)[2])
  } else {
    mu = rnorm(m)
  }
  lambda = rgamma(K,a0,b0)
  beta = rnorm(zm)
  if (ar1priorBETA==1){
    ar1 = ifelse(true$structured, rbeta(1,ar11,ar12),0.5)
  } else {
    ar1 = ifelse(true$structured, runif(1,-1,1),0.5)
  }
  alpha = rnorm(m,0,0.01)
  if (true$structured) {
    r12 = rnorm(choose(M,2))
    Sigma = ar1_cor(J,ar1) %x% Corrmat(M,r12)$Gamma
    Gammainv = matrix(NA, m+1,m+1) ; Gammainv[1,1] = crossprod(alpha, solve(Sigma,alpha)) + 1 ; Gammainv[-1,-1] = solve(Sigma) ; Gammainv[-1,1] = Gammainv[1,-1] = -solve(Sigma,alpha)
    Gamma = solve(Gammainv)
  } else {
    r12   = rnorm(choose(m+1,2))
    Gamma = Corrmat(m+1,r12)$Gamma
  }
  
}
if (!true$structured) Sigma = diag(m)
Gammamm = Gamma[2:(m+1),2:(m+1)] ; Gammainv = solve(Gamma) ; Gammamminv = solve(Gammamm)
alpha = alphat
inits = list(sig2=rep(sig2,M), mu=mu, Sigma=Sigma, Gamma=Gamma, r12=r12, lambda=lambda, beta=beta, ar1=ar1, alpha=alpha, Gammastr_ar1=1*Gammastr_ar1)

### SLICE SAMPLING PARAMETER ###
wlong = ifelse(sd(Y[,1], na.rm=TRUE) > 1,1,sd(Y[,1], na.rm=TRUE))

structured = 1*structured ; MAR = 1*MAR ; ar1priorBETA = 1*ar1priorBETA


MM = 10000 ; burnin = 1000 ; thin = 5 ; saving = 1000 ; MHmcmc=500 ; accepttol = 0.35
priors = dd$priors ; slice_param = list(m=1,w=wlong)
dim = list(J=M,M=J) ; structured = true$structured

start = Sys.time()

postsampling = MCMC_pwe(priors=dd$priors, Y = dat[,grep('Y',colnames(dat))], X=as.matrix(dat[,grep('X',colnames(dat))]), Z=as.matrix(dat[,grep('Z',colnames(dat))]), v = dat[,"v"], obsind = dat[,"obsind"], Tst=dat[,"Tst"], delta=dat[,"delta"], dat=dat, W=W, inits=inits, timeint=timeint, dim=list(J=M,M=J), slice_param=list(m=1,w=wlong), MAR=MAR, structured=true$structured, IMidx=IMidx, M=MM, burnin=burnin,  MHmcmc=MHmcmc, thin=thin, saving=saving, accepttol=accepttol, filename=paste0(node.idx,"Simdat",nnn,"-chain",nchain,"PWEcpp_str_",true$structured))

end = Sys.time()

colnames(postsampling$dat) = colnames(dat)
postsampling$time = end-start ; postsampling$MLE = MLE.estimates
filename=paste0(node.idx,"Simdat",nnn,"-chain",nchain,"PWEcpp_str_",true$structured,filename_ar1prior,"MAR_",MAR)
saveRDS(postsampling, file=paste0("./results/",format(Sys.time(),"%y%m%d_%H%M%S"),"-",filename,".RData"))

