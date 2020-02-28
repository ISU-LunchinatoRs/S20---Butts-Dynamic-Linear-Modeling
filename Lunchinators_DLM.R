# Lunchinators -- Online Dynamic Linear Modeling 

rm(list = ls())
# Application of the Online Dynamic Linear Model 

setwd("C:/Users/Owner/Desktop/R/Lunchinators/DLM analysis") # set wd to wherever you have the files

## Peter Lake Example -------------------------------------------------

# Function to fit an AR(p) model as a Dynamic Linear Model (DLM) by Maximum Likelihood
# Copyright 2017 by Stephen R. Carpenter

# The DLM, a form of dynamic linear regression, is defined by

# Observation equation is
# Y_t = F_t'*theta_t + eta_t where
# Y_t is the prediction
# F_t is a vector of predictors at the beginning of the time step
# theta_t is the parameter vector
# eta_t is an individual observation error

# System equation is:
# theta_t = theta_t-1 + omega_t
# where theta is defined above and omega_t is an individual process error

# This model converts directly into MARSS format as follows: 
#  x(t) = B(t)*x(t-1) + w(t)
#  y(t) = Z(t)*x(t) + v(t)
# where x(t) = theta_t, B(t) = G(t) = I, Z(t) = F_t, Q=Q, R=r, and
#   u=0, C=0, c=0 a=0, D=0, d=0
# For more information about MARSS see the user guide; the version I used was
#
# Holmes EE, Ward EJ, Scheuerell MD, 2014, Analysis of Multivariate Time Series Using the
# MARSS package. V 3.9. Northwest Fisheries Science Center, NOAA, Seattle, WA, USA
#
# Function call: DLM.MLE(nl,nobs,timevec,X,title) where
#  nl is the number of AR lags, nobs is the number of observations, timevec is the time steps,
#  X is the time series to be analyzed, and title is the title for the plots

# Outputs are a list containing:
#  1 = matrix containing: time step, Y, yhat (one-step prediction); 
#       Dimension is (nobs-nl)x3 where nobs is number of observations and nl is number of lags
#  2 = (nobs-nl)x4 matrix containing eigenvalue, sd of eigenvalue, eigenvalue + sd, eigenvalue - sd
#  3 = (nl+1)x(nobs-nl) matrix of AR parameter estimates; col 1 is intercept, col 2 is AR(1) coef, etc.
#  4 = (nl+1)x(nobs-nl) matrix of AR parameter standard deviations
#  5 = the full MARSS structure returned by MARSS


# install.packages('MARSS') #Uncomment if need to install package 
library(MARSS)

# MLE fit of DLM using MARSS: input nl, nobs, X, timevec *************************************************
DLM.MLE = function(nl,nobs,timevec,X,title)  {
  
  if(nl>8) {
    nl=8
    print('Number of lags nl must be 8 or less; nl set to 8',quote=F)
    print('High order AR fits are likely to have convergence problems',quote=F)
  }
  
  # number of parameters
  p = nl+1 # number of parameters counting intercept
  
  # AR variates
  X.design = matrix(1,nr=(nl+1),nc=(nobs-nl)) # matrix to hold predictors
  for(i in 1:nl) {
    X.design[(i+1),] = X[i:(nobs-nl+i-1)]
  }
  Y = matrix(X[(1+nl):nobs],nr=1,nc=(nobs-nl)) # response
  
  # LS regression
  invXX = solve(X.design%*%t(X.design))
  lm.par = invXX%*%X.design%*%t(Y)
  # Force initials inside unit circle
  lm.inits=lm.par
  lm.inits[2:p] = ifelse(lm.par[2:p]^2 < 1,lm.par[2:p],0.9)
  
  # process equation
  B = diag(p)
  U = matrix(0,nr=p,nc=1) 
  Q = matrix(list(0),p,p)
  vec.of.qs = c('q1','q2','q3','q4','q5','q6','q7','q8','q9')
  for(i in 1:p) {
    Q[i,i] = vec.of.qs[i]
  }
  
  # Observation equation
  Z = array(NA,c(1,p,(nobs-nl))) # first dimension is time, second is parameters, third is time steps
  for(i in 1:p) {
    Z[1,i,] = X.design[i,]
  }
  A = matrix(0)
  R = matrix('r')
  
  # initials
  inits.list = list(x0=matrix(lm.inits, nrow=p))
  
  # model list
  mod.list = list(B=B, U=U, Q=Q, Z=Z, A=A, R=R)
  
  # run DLM
  dlm1 = MARSS(Y, inits=inits.list, model=mod.list,control=list(maxit=5000,trace=1))
  
  # extract state estimates
  B.ests = dlm1$states
  # s.e. of state estimates
  B.se = dlm1$states.se
  
  # one-step predictions
  Yhat=rep(0,(nobs-nl))
  for(i in 1:(nobs-nl)){
    Yhat[i] = X.design[,i]%*%B.ests[,i]
  }
  
  epsilon = Yhat-Y # prediction errors
  
  # compute eigenvalues
  lamda = rep(0,(nobs-nl))
  for(i in 1:(nobs-nl)) {
    armat = matrix(0,nr=nl,nc=nl)
    subdiag = rep(1,(nl-1))
    armat[row(armat) == col(armat)+1] = subdiag
    armat[1,] = B.ests[2:p,i]
    eigvals = eigen(armat,only.values=T)
    lamda[i] = max(Mod(eigvals$values))
  }
  
  # Bootstrap eigenvalue errors
  NT = nobs-nl  # number of time steps
  nboot = 100 # number of bootstrap iterations
  lamda.se = rep(0,NT) # vector to hold time series of eigenvalue s.e.
  lamda.temp = rep(0,nboot)  # vector to hold temporary sample of eigenvalues
  # make a matrix with subdiagonal 1
  armat = matrix(0,nr=nl,nc=nl)
  subdiag = rep(1,(nl-1))
  armat[row(armat) == col(armat)+1] = subdiag
  for(i in 1:NT) {  # loop over time steps
    # make a matrix of normal random numbers with nboot rows and nl columns, each column a sample for an AR coef
    ARnorm = mapply(function(mu,sigma){rnorm(n=nboot,mean=mu,sd=sigma)},mu=B.ests[2:p,i],sigma=B.se[2:p,i])
    for(j in 1:nboot)  {  # loop over bootstrap iterations
      armat[1,] = ARnorm[j,]
      eigvals = eigen(armat,only.values=T)
      lamda.temp[j] = max(Mod(eigvals$values))
    }  # end loop over bootstrap iterations
    lamda.se[i] = sd(lamda.temp)
  }  # end loop over timesteps
  
  # Time vector for plotting
  T.ar = timevec[1:(nobs-nl)]
  
  # Plot 1: general results
  windows(width=6,height=12)
  par(mfrow=c(4,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
  plot(T.ar,Yhat,type='l',lwd=3,col='deepskyblue',xlab='DOY',ylab='Y and Yhat',main=title)
  points(T.ar,Y,type='p',pch=19,col='darkblue')
  plot(T.ar,epsilon,type='b',pch=19,lwd=2,col='red',xlab='DOY',ylab='one-step error')
  plot(T.ar,B.ests[1,],type='b',pch=19,lwd=2,col='forestgreen',xlab='DOY',ylab='intercept')
  eigrange = range(lamda,1)
  plot(T.ar,lamda,type='b',pch=19,lwd=2,col='blue',xlab='DOY',ylab='eigenvalue',ylim=eigrange)
  abline(h=1,lty=2)
  
  # Plot 2: plot the ar coefficients
  windows(width=6,height=12)
  par(mfrow=c((nl),1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
  for(i in 2:p) {  # plot AR coefs
    plot(T.ar,B.ests[i,],type='l',lwd=2,col='magenta',
         ylab=bquote('AR'~ .(i-1) ~'coef'),
         xlab='DOY')
  }
  
  # Plot 3: eigenvalue plus error
  lamda.plus = lamda+lamda.se
  lamda.minus = lamda-lamda.se
  windows()
  par(mfrow=c(1,1),mar=c(4, 4.3, 2, 2) + 0.1, cex.axis=1.6,cex.lab=1.6)
  yrange = range(lamda.plus,lamda.minus,0,1)
  plot(T.ar,lamda,type='l',lwd=2,col='blue',ylim=yrange,
       xlab='DOY',ylab='Eigenvalue +/- SE',main=title)
  points(T.ar,lamda.plus,type='l',lwd=2,lty=2,col='deepskyblue')
  points(T.ar,lamda.minus,type='l',lwd=2,lty=2,col='deepskyblue')
  abline(h=1,lty=2)
  
  # Construct output list
  Yyhat = matrix(c(T.ar,Y,Yhat),nr=(nobs-nl),nc=3)
  LamdaMat = matrix(c(lamda,lamda.se,lamda.plus,lamda.minus),nr=(nobs-nl),nc=4)
  
  outlist = list(Yyhat,LamdaMat,B.ests,B.se,dlm1)
  return(outlist)
  
}  # End function for fitting DLM with MLE **********************************************************************



# Try time-varying AR on Cascade data using ML DLM method with MARSS
# SRC 2017-07-01
# (c) Stephen R. Carpenter 2017-07-01


# Load Peter Lake demo data set
# save(x.full,T.full,file='PeterL_2015_manual_chl.Rdata')
load(file='PeterL_2015_BGA_HYLB.Rdata')

# log transform? Optional. Comment out to use natural unit
x.full = log10(x.full)
title = c('Peter 2015 log10 Phycocyanin (HYLB), ML DLM AR(p)')
nobs = length(x.full)

# Detrend and standardize (optional)
linmod = lm(x.full ~ T.full)
err.std = linmod$residuals/sd(linmod$residuals)
x.full = err.std  # comment out this line to NOT detrend

# set up for function call
nl = 2 # number of lags = AR order
X = x.full
timevec = T.full

# Compute dlm
dlmfit = DLM.MLE(nl,nobs,timevec,X,title)
Yyhat = dlmfit[[1]]
EigenVals = dlmfit[[2]]
B.est = dlmfit[[3]]
B.se = dlmfit[[4]]
dlm.struc = dlmfit[[5]]

# Outputs are a list containing:
#  1 = matrix containing: time step, Y, yhat (one-step prediction); 
#       Dimension is (nobs-nl)x3 where nobs is number of observations and nl is number of lags
#  2 = (nobs-nl)x4 matrix containing eigenvalue, sd of eigenvalue, eigenvalue + sd, eigenvalue - sd
#  3 = (nl+1)x(nobs-nl) matrix of AR parameter estimates; col 1 is intercept, col 2 is AR(1) coef, etc.
#  4 = (nl+1)x(nobs-nl) matrix of AR parameter standard deviations
#  5 = the full MARSS structure returned by MARSS

## Hort Farm example --------------------------------

rm(list = ls())

## MODIFIED ONLINE DYNAMIC LINEAR MODEL (DLM) ESTIMATION 

DLM <- function(delta,n.gamma,d.gamma,mvec,Cpar,Yvec,Fmat) {
  
  # Online algorithm for Dynamic linear regression
  # Copyright 2016 by Stephen R. Carpenter
  
  # Description and definitions:
  
  # Observation equation is
  # Y_t = F_t'*theta_t + eta_t where
  # Y_t is the prediction
  # F_t is a vector of predictors at the beginning of the time step
  # theta_t is the parameter vector
  # eta_t is an individual observation error
  
  # System equation is:
  # theta_t = theta_t-1 + omega_t
  # where theta is defined above and omega_t is an individual process error
  
  # Inputs to the function are:
  # delta, the discount factor
  # n.gamma, the initial number of observations (usually 1)
  # d.gamma, the initial shape parameter for prediction errors
  #  (prior estimate of prediction variance = d.gamma / n.gamma)
  # mvec, the initial guess of regression coefficients
  # Cpar, the initial guess of the covariance matrix of regression coefficients
  # Yvec, the vector of the observed response variate
  # Fmat, the matrix of predictors
  
  # Outputs are:
  # predix, the one-step-ahead predictions of the response variate
  # varpredix, the prediction variance at start of time step before error is measured
  # pars, the updated parameter estimates using the most recent prediction error
  # parvar, the variances of the parameters
  # Svec, the update (after error is measured within a time step) of varpredix
  
  # Updating follows the equations on p. 176-179 of Carpenter 2003,
  # Regime Shifts in Lake Ecosystems: Pattern and Variation
  
  # Determine constants
  npar <- length(mvec)
  Nobs <- length(Yvec)
  S0 <- d.gamma/n.gamma
  
  # Set up vectors to hold results
  predix <- rep(0,Nobs)
  varpredix <- rep(0,Nobs)
  Svec = rep(0,Nobs)
  pars <- matrix(0,nrow=Nobs,ncol=npar)
  parvar = matrix(0,nrow=Nobs,ncol=npar)
  
  for(i in 1:Nobs)  {  #Start DLM loop
    # Generate predictions
    Fvec <- Fmat[i,] # vector of predictors
    predix[i] <- sum(Fvec*mvec)
    # Compute error and update estimates
    error <- Yvec[i]-predix[i]
    Rmat <- Cpar/delta
    varpredix[i] <- (t(Fvec) %*% Rmat %*% Fvec) + S0
    n.gamma <- (delta*n.gamma)+1
    d.gamma <- (delta*d.gamma)+(S0*error*error/varpredix[i])
    S1 <- d.gamma/n.gamma
    Svec[i] = S1  # save updated variance
    Avec <- (Rmat %*% Fvec)/varpredix[i]
    mvec <- mvec + (Avec*error)
    pars[i,] <- mvec
    Cpar <- (S1/S0)*(Rmat - (Avec %*% t(Avec))*varpredix[i])
    # Disallow negative variances on the diagonal
    for(idiag in 1:npar) {
      Cpar[idiag,idiag] <- max(0,Cpar[idiag,idiag])
    }
    parvar[i,] = diag(Cpar)
    S0 <- S1 # roll over S
  } # End DLM loop
  
  DLM.out <- list(predix,varpredix,pars,parvar,Svec)
  return(DLM.out)
} # END DLM FUNCTION

# ONLINE SHELL 

ODLMAR = function(nl,delta,x.full,T.full,title) {
  # Compute online DLM for AR models
  
  # Copyright 2016 by Stephen R. Carpenter
  
  # This function is a shell for the DLM() function
  
  # Inputs are:
  #  nl = number of lags in AR model
  #  delta = discount factor 0<delta<1; reasonable values are 0.9 to 0.99
  #     Rule of Thumb: df for each point estimate = 1/(1-delta)
  #  x.full is the time series to be analyzed
  #  T.full is the corresponding time steps
  #  title is a title for the plots
  
  # Outputs are a list containing:
  #  1 = matrix containing: time step, Y, yhat (one-step prediction), updated prediction variance
  #     Dimension is (nobs-nl)x4 where nobs is number of observations and nl is number of lags
  #  2 = (nobs-nl)x4 matrix containing eigenvalue, sd of eigenvalue, eigenvalue + sd, eigenvalue - sd
  #  3 = (nl+1)x(nobs-nl) matrix of AR parameter estimates; col 1 is intercept, col 2 is AR(1) coef, etc.
  #  4 = (nl+1)x(nobs-nl) matrix of AR parameter standard deviations
  
  # choose AR order
  p = nl+1 # allow for intercept
  
  # Number of observations
  nobs = length(x.full)
  
  # AR variates
  X.design = matrix(1,nr=(nl+1),nc=(nobs-nl)) # matrix to hold predictors
  for(i in 1:nl) {
    X.design[(i+1),] = x.full[i:(nobs-nl+i-1)]
  }
  Y = matrix(x.full[(1+nl):nobs],nr=1,nc=(nobs-nl)) # response
  
  # LS regression for initials
  invXX = solve(X.design%*%t(X.design))
  lm.par = invXX%*%X.design%*%t(Y)
  # Force initial parameters inside unit circle
  lm.inits=lm.par
  lm.inits[2:p] = ifelse(lm.par[2:p]^2 < 1,lm.par[2:p],0.9)
  # Other useful regression statistics
  lm.yhat = t(X.design)%*%lm.par
  lm.err = t(Y)-lm.yhat
  verr = var(lm.err)
  covpar = invXX*verr[1,1]
  
  # Other parameters for DLM
  n.gam = 1  # initial df
  d.gam = verr[1,1]  # based on the gamma distribution mean var = df*scale
  
  # Run DLM
  DLMrun = DLM(delta,n.gam,d.gam,lm.par,covpar,Y,t(X.design))
  yhat = DLMrun[[1]]
  vyhat = DLMrun[[2]]
  B.ests = t(DLMrun[[3]])
  B.sd = t( sqrt(DLMrun[[4]]) ) # parameter SD
  vupdate = DLMrun[[5]] #updated variance
  
  # Compute adjusted R^2
  Ry = cor(t(Y),yhat)
  R2 = Ry^2
  R2adj = 1 - ( (1-R2[1,1])*(nobs-1)/(nobs-p-1) )
  print('',quote=F)
  print(c('DLM R^2 = ',round(R2,3),' DLM adj R^2 = ',round(R2adj,3)),quote=F)
  # Compute AIC
  err1 = yhat-Y # errors
  sd1 = sqrt(vupdate)  # error variance estimate
  LL = dnorm(err1, mean=0, sd=sd1, log=T) # log likelihoods
  aic = 2*p - 2*sum(LL)
  print(c('AIC = ',round(aic,2)),quote=F)
  
  # compute eigenvalues
  lamda = rep(0,(nobs-nl))
  for(i in 1:(nobs-nl)) {
    armat = matrix(0,nr=nl,nc=nl)
    subdiag = rep(1,(nl-1))
    armat[row(armat) == col(armat)+1] = subdiag
    armat[1,] = B.ests[2:p,i]
    eigvals = eigen(armat,only.values=T)
    lamda[i] = max(Mod(eigvals$values))
  }
  
  # Bootstrap eigenvalue errors
  NT = nobs-nl  # number of time steps
  nboot = 100 # number of bootstrap iterations
  lamda.sd = rep(0,NT) # vector to hold time series of eigenvalue s.e.
  lamda.temp = rep(0,nboot)  # vector to hold temporary sample of eigenvalues
  # make a matrix with subdiagonal 1
  armat = matrix(0,nr=nl,nc=nl)
  subdiag = rep(1,(nl-1))
  armat[row(armat) == col(armat)+1] = subdiag
  for(i in 1:NT) {  # loop over time steps
    # make a matrix of normal random numbers with nboot rows and nl columns, each column a sample for an AR coef
    ARnorm = mapply(function(mu,sigma){rnorm(n=nboot,mean=mu,sd=sigma)},mu=B.ests[2:p,i],sigma=B.sd[2:p,i])
    for(j in 1:nboot)  {  # loop over bootstrap iterations
      armat[1,] = ARnorm[j,]
      eigvals = eigen(armat,only.values=T)
      lamda.temp[j] = max(Mod(eigvals$values))
    }  # end loop over bootstrap iterations
    lamda.sd[i] = sd(lamda.temp)
  }  # end loop over timesteps
  
  # Plot results
  T.ar = T.full[1:(nobs-nl)]
  
  # Plot 3: eigenvalue plus error
  lamda.plus = lamda+lamda.sd
  lamda.minus = lamda-lamda.sd
  yrange = range(lamda.plus,lamda.minus,0,1)
  plot(T.ar,lamda,type='l',lwd=2,col=color,ylim=c(0,4),
       xlab='',ylab='',main=title, cex.axis=1.5)
  polygon(c(T.ar, rev(T.ar)), c(lamda.plus, rev(lamda.minus)), col=color1, border=NA)
  points(T.ar,lamda,type='l',lwd=4,col=color)
  # points(T.ar,lamda.plus,type='l',lwd=2,lty=2,col=color1)
  # points(T.ar,lamda.minus,type='l',lwd=2,lty=2,col=color)
  abline(h=1,lty=2)
  abline(v=170,lty=1,lwd=2)
  #Lines for Heat Wave
  #abline(v=190,lty=2, lwd=2)
  #abline(v=201,lty=2, lwd=2)
  
  
  # Create output list and return it
  Yyhat = matrix(c(T.ar,Y,yhat,vupdate),nr=(nobs-nl),nc=4)
  LamdaMat = matrix(c(lamda,lamda.sd,lamda.plus,lamda.minus),nr=(nobs-nl),nc=4)
  
  outlist = list(Yyhat,LamdaMat,B.ests,B.sd)
  return(outlist)
  
  
}  # End online DLM function

# end functions 

# Hort Farm Data 
  # Make sure the ODLMAR and DLM function is run

# High Planktivore Pond 

bmb19 = read.csv("DLM_DAT.csv")
x.full = as.vector(bmb19[bmb19$POND=='F' & bmb19$FISH=="post", "CHL"])
T.full = as.vector(bmb19[bmb19$POND=="F" & bmb19$FISH=="post", "DOY"])

plot(T.full, x.full)

# log transform? Optional. Comment out to use natural unit
# x.full = log10(x.full)
title = c('Pond F - 2019, Online DLM AR(1)') #Title for the plot
color= "#5ab4ac" # this is the color for the eigenvalue line
color1="gray88" #this is the color for the error polygon 
nobs = length(x.full)

# START PROTOTYPE SHELL
# USER MUST INPUT: nl; delta; x.full; T.full; title

nl = 1 # number of lags, for Chlorophyll we decided AR(1) based on AIC
delta = 0.9 # 0<delta<1; see advice in functions

ODL.out = ODLMAR(nl,delta,x.full,T.full,title)

Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]
B.sd = ODL.out[[4]]

#COLORS 
# "#d8b365" = control
# "#5ab4ac" = pulse 

# Reference High Planktivore Pond 
bmb19 = read.csv("DLM_DAT.csv")
x.full = as.vector(bmb19[bmb19$POND=='B' & bmb19$FISH=="post", "CHL"])
T.full = as.vector(bmb19[bmb19$POND=="B" & bmb19$FISH=="post", "DOY"])

plot(T.full, x.full)

# log transform? Optional. Comment out to use natural unit
# x.full = log10(x.full)
title = c('Pond B - 2019, Online DLM AR(1)') #Title for the plot
color= "#d8b365" # this is the color for the eigenvalue line
color1="gray88" #this is the color for the error polygon 
nobs = length(x.full)

# START PROTOTYPE SHELL
# USER MUST INPUT: nl; delta; x.full; T.full; title

nl = 1 # number of lags, for Chlorophyll we decided AR(1) based on AIC
delta = 0.9 # 0<delta<1; see advice in functions

ODL.out = ODLMAR(nl,delta,x.full,T.full,title)

Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]
B.sd = ODL.out[[4]]

# Low Planktivore Pond 
bmb19 = read.csv("DLM_DAT.csv")
x.full = as.vector(bmb19[bmb19$POND=='E' & bmb19$FISH=="post", "CHL"])
T.full = as.vector(bmb19[bmb19$POND=="E" & bmb19$FISH=="post", "DOY"])

plot(T.full, x.full)

# log transform? Optional. Comment out to use natural unit
# x.full = log10(x.full)
title = c('Pond E - 2019, Online DLM AR(1)') #Title for the plot
color= "#5ab4ac" # this is the color for the eigenvalue line
color1="gray88" #this is the color for the error polygon 
nobs = length(x.full)

# START PROTOTYPE SHELL
# USER MUST INPUT: nl; delta; x.full; T.full; title

nl = 1 # number of lags, for Chlorophyll we decided AR(1) based on AIC
delta = 0.9 # 0<delta<1; see advice in functions

ODL.out = ODLMAR(nl,delta,x.full,T.full,title)

Yyhat = ODL.out[[1]]
EigenVals = ODL.out[[2]]
B.ests = ODL.out[[3]]
B.sd = ODL.out[[4]]




