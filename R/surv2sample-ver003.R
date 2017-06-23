#' @name surv2sampleComp-package
#' @aliases  surv2sampleComp-package
#' @docType  package
#' @title Inference for Model-Free Between-Group Parameters For Censored Survival Data
#' @description
#' Performs inference of several model-free group contrast measures, which include difference/ratio of cumulative incidence rates
#' at given time points, quantiles, and restricted mean survival times (RMST).
#' Two kinds of covariate adjustment procedures (i.e., regression and augmentation) for inference of the metrics based on RMST are also included.
#' @author Lu Tian, Hajime Uno, Miki Horiguchi
#'
#' Maintainer: Miki Horiguchi <horiguchimiki@gmail.com>
#'
#' @references
#' Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.
#'
#' Zhao L, Tian L, Uno H, Solomon S, Pfeffer M, Schindler J, Wei LJ. Utilizing the integrated difference of two survival functions to quantify the treatment contrast for designing, monitoring, and analyzing a comparative clinical study. Clinical Trials 2012, 9, 570-577.
#'
#' @keywords
#' survival
#' @seealso flexsurv plotrix survival
#' @import flexsurv plotrix survival
#' @importFrom graphics plot
#' @importFrom stats glm lm integrate pchisq pnorm qnorm quantile rexp sd
#' @importFrom utils data
NULL


#' @name pbc.sample
#' @aliases  pbc.sample
#' @title Edit pbc data to run sample code
#' @description Edit pbc data in survival package and make it ready to run the sample code in this manual.
#' @usage pbc.sample()
#' @seealso \code{pbc} in survival package
NULL


CompCase=function(mydata){
  sum(is.na(apply(mydata,1,mean)))
  mydata<-mydata[!is.na(apply(mydata,1,mean)),]
}
NULL

#' @export
#######################################
# pbc.sample
######################################
pbc.sample=function(){
  Z=list()
  D=CompCase(survival::pbc[,c(2:4,10:14)]);
  Z$time=D$time/365.25;
  Z$status=as.numeric(D$status==2);
  Z$group=as.numeric(D$trt==2);
  Z$covariates=as.matrix(D[,4:8]);
  Z
}
NULL


#' @name plot.surv2sample
#' @aliases plot.surv2sample
#' @title Plot method for surv2sample objects
#' @description Creates plots from a surv2sample object.
#' @param x surv2sample object
#' @param measure The type of measure used for the plot. When default(=NULL), plot.survfit() is called and KM plots are given.
#' When "relative time" is specified, a plot of relative percentiles with corresponidng 0.95 confidence intervals is generatead.
#' @param baseline Indicates the baseline group, 0/1. Default is 0.
#' @param ... For further method
#' @seealso \code{plotCI} in plotrix package
#' @export
######################################
# plot.surv2sample
######################################
plot.surv2sample=function(x, measure=NULL,baseline=0,...){

  if(is.null(measure)){
    y=x$survfit ; class(y)="survfit" ; plot(y,...)
  }

  if(measure=="relative percentile"){
    xx=x$quanprobs
    if(baseline==0){y=x$contrast.ratio01}else{y=x$contrast.ratio10}
    yy=y[(nrow(y)-length(xx)+1):nrow(y),]
    plotrix::plotCI(xx, yy[,1], uiw=yy[,3], liw=yy[,2], xlab="percent", ylab="relative time (95%CI)",...)
  }

}
NULL


#' @name rmstreg
#' @aliases rmstreg
#' @title Adjusted difference/ratio of restricted mean survival times
#' @description Compares restricted mean survival time between two groups, adjusting for imbalance of baseline factors via a regression model.
#' @usage rmstreg(y, delta, x, arm, tau, type="difference", conf.int=0.95)
#' @param y The follow-up time.
#' @param delta The censoring indicator, 1=event, and 0=censoring.
#' @param x The covariate matrix. The first colomn of this matrix should be the group indicator, arm (below).
#' @param arm The group indicator, 1/0.
#' @param tau The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time.
#' @param type The type of the between-group contrast measure: "difference"(default), "ratio" or "lossratio".
#' @param conf.int The level for computation of the confidence intervals. The default is 0.95.
#' @author Lu Tian
#' @references
#' Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.
#' @examples
#' D=pbc.sample()
#' x=cbind(D$group, D$covariates)
#' rmstreg(D$time, D$status, x, D$group, tau=8, type="difference")
#' @export
######################################################################
#included in ainternalmedicine-ver001.R (ver 1.0-1)
rmstreg=function(y, delta, x, arm, tau, type="difference", conf.int=0.95)
{
  label.low = paste0("lower ", conf.int)
  label.upp = paste0("upper ", conf.int)

  if(type!="difference" && type!="ratio" && type!="lossratio")
    print("Type must be difference, ratio or lossratio.")

  if(type=="difference" || type=="ratio" || type=="lossratio"){

    n=length(y)
    x=cbind(1, x)
    p=length(x[1,])

    y0=pmin(y, tau)
    d0=delta
    d0[y0==tau]=1

    d10=d0[arm==1]
    d00=d0[arm==0]
    y10=y0[arm==1]
    y00=y0[arm==0]
    x1=x[arm==1,]
    x0=x[arm==0,]
    n1=length(d10)
    n0=length(d00)


    id1=order(y10)
    y10=y10[id1]
    d10=d10[id1]
    x1=x1[id1,]

    id0=order(y00)
    y00=y00[id0]
    d00=d00[id0]
    x0=x0[id0,]

    fitc1=survfit(Surv(y10, 1-d10)~1)
    fitc0=survfit(Surv(y00, 1-d00)~1)

    weights1=d10/rep(fitc1$surv, table(y10))
    weights0=d00/rep(fitc0$surv, table(y00))

    weights=c(weights1, weights0)

    if(type=="difference")
    {fitt=lm(c(y10,y00)~rbind(x1, x0)-1, weights=weights)
    beta0=fitt$coef

    error1=y10-as.vector(x1%*%beta0)
    score1=x1*weights1*error1

    error0=y00-as.vector(x0%*%beta0)
    score0=x0*weights0*error0
    }

    if(type=="ratio")
    {fitt=glm(c(y10,y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=y10-exp(as.vector(x1%*%beta0))
    score1=x1*weights1*error1

    error0=y00-exp(as.vector(x0%*%beta0))
    score0=x0*weights0*error0
    }


    if(type=="lossratio")
    {fitt=glm(c(tau-y10,tau-y00)~rbind(x1, x0)-1, family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=tau-y10-exp(as.vector(x1%*%beta0))
    score1=x1*weights1*error1

    error0=tau-y00-exp(as.vector(x0%*%beta0))
    score0=x0*weights0*error0
    }



    kappa.arm1=matrix(0, n1, p)
    for(i in 1:n1)
    {kappa1=score1[i,]

    kappa2=apply(score1[y10>=y10[i],,drop=F], 2, sum)*(1-d10[i])/sum(y10>=y10[i])

    kappa3=rep(0, p)
    for(k in 1:n1)
    { if(y10[k]<=y10[i])
      kappa3=kappa3+apply(score1[y10>=y10[k],,drop=F], 2, sum)*(1-d10[k])/(sum(y10>=y10[k]))^2
    }

    kappa.arm1[i,]=kappa1+kappa2-kappa3
    }


    kappa.arm0=matrix(0, n0, p)
    for(i in 1:n0)
    {kappa1=score0[i,]

    kappa2=apply(score0[y00>=y00[i],,drop=F], 2, sum)*(1-d00[i])/sum(y00>=y00[i])

    kappa3=rep(0, p)
    for(k in 1:n0)
    {if(y00[k]<=y00[i])
      kappa3=kappa3+apply(score0[y00>=y00[k],,drop=F], 2, sum)*(1-d00[k])/(sum(y00>=y00[k]))^2
    }

    kappa.arm0[i,]=kappa1+kappa2-kappa3
    }


    if(type=="difference")
    {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
    A=t(x)%*%x
    varbeta=solve(A)%*%gamma%*%solve(A)
    }

    if(type=="ratio" || type=="lossratio")
    {gamma=t(kappa.arm1)%*%kappa.arm1+t(kappa.arm0)%*%kappa.arm0
    A=t(x*exp(as.vector(x%*%beta0)))%*%x
    varbeta=solve(A)%*%gamma%*%solve(A)
    }


    if(type=="difference")
    {beta0=beta0
    se0=sqrt(diag(varbeta))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    cilow=beta0-se0*abs(qnorm((1-conf.int)/2))
    cihigh=beta0+se0*abs(qnorm((1-conf.int)/2))
    result=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, label.low=cilow, label.upp=cihigh)
    colnames(result)=c("coef", "se(coef)", "z", "p", label.low, label.upp)
    }

    if(type=="ratio" || type=="lossratio")
    {beta0=beta0
    se0=sqrt(diag(varbeta))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    r0=exp(beta0)
    cilow=exp(beta0-se0*abs(qnorm((1-conf.int)/2)))
    cihigh=exp(beta0+se0*abs(qnorm((1-conf.int)/2)))
    result=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "exp(coef)"=exp(beta0), label.low=cilow, label.upp=cihigh)
    colnames(result)=c("coef", "se(coef)", "z", "p", "exp(coef)", label.low, label.upp)
    }

    if(p==2)
      rownames(result)=c("intercept", "x")

    if(p>2)
      rownames(result)=c("intercept", colnames(x[,-1]))

    return(result)
  }
}
NULL


#' @name rmstaug
#' @aliases rmstaug
#' @title Adjusted difference/ratio of restricted mean survival times
#' @description Compares restricted mean survival time between two groups, adjusting for imbalance of baseline factors.
#' @usage rmstaug(y, delta, x, arm, tau, type="difference", conf.int=0.95)
#' @param y The follow-up time.
#' @param delta The censoring indicator, 1=event, and 0=censoring.
#' @param x The covariate matrix. The group indicator, arm (below) should not be included in this matrix.
#' @param arm The group indicator, 1/0.
#' @param tau The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time.
#' @param type The type of the between-group contrast measure: "difference"(default), "ratio" or "lossratio".
#' @param conf.int The level for computation of the confidence intervals. The default is 0.95.
#' @author Lu Tian
#' @references
#' Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.
#' @examples
#' D=pbc.sample()
#' rmstaug(D$time, D$status, D$covariates, D$group, tau=8, type="difference")
#' @export
##################################################################################
#included in ainternalmedicine-ver001.R (ver 1.0-1)
rmstaug=function(y, delta, x, arm, tau, type="difference", conf.int=0.95)
{
  label.low = paste0("lower ", conf.int)
  label.upp = paste0("upper ", conf.int)

  if(type!="difference" && type!="ratio" && type!="lossratio")
    print("Type must be difference, ratio or lossratio.")

  if(type=="difference" || type=="ratio" || type=="lossratio"){

    n=length(y)
    x=as.matrix(x)
    p=length(x[1,])
    pi=mean(arm)


    y0=pmin(y, tau)
    d0=delta
    d0[y0==tau]=1

    d10=d0[arm==1]
    d00=d0[arm==0]
    y10=y0[arm==1]
    y00=y0[arm==0]
    x1=x[arm==1,,drop=F]
    x0=x[arm==0,,drop=F]
    n1=length(d10)
    n0=length(d00)


    id1=order(y10)
    y10=y10[id1]
    d10=d10[id1]
    x1=x1[id1,,drop=F]

    id0=order(y00)
    y00=y00[id0]
    d00=d00[id0]
    x0=x0[id0,,drop=F]

    fitc1=survfit(Surv(y10, 1-d10)~1)
    fitc0=survfit(Surv(y00, 1-d00)~1)

    weights1=d10/rep(fitc1$surv, table(y10))
    weights0=d00/rep(fitc0$surv, table(y00))

    weights=c(weights1, weights0)

    if(type=="difference")
    {fitt=lm(c(y10,y00)~rep(c(1, 0), c(n1, n0)), weights=weights)
    beta0=fitt$coef

    error1=y10-beta0[1]-beta0[2]
    score1=cbind(1, rep(1, n1))*weights1*error1

    error0=y00-beta0[1]
    score0=cbind(1, rep(0, n0))*weights0*error0
    }

    if(type=="ratio")
    {fitt=glm(c(y10,y00)~rep(c(1, 0), c(n1, n0)), family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=y10-exp(beta0[1]+beta0[2])
    score1=cbind(1, rep(1, n1))*weights1*error1

    error0=y00-exp(beta0[1])
    score0=cbind(1, rep(0, n0))*weights0*error0
    }


    if(type=="lossratio")
    {fitt=glm(c(tau-y10,tau-y00)~rep(c(1, 0), c(n1, n0)), family="quasipoisson", weights=weights)
    beta0=fitt$coef

    error1=tau-y10-exp(beta0[1]+beta0[2])
    score1=cbind(1, rep(1, n1))*weights1*error1

    error0=tau-y00-exp(beta0[1])
    score0=cbind(1, rep(0, n0))*weights0*error0
    }



    kappa.arm1=matrix(0, n1, 2)
    for(i in 1:n1)
    {kappa1=score1[i,]

    kappa2=apply(score1[y10>=y10[i],,drop=F], 2, sum)*(1-d10[i])/sum(y10>=y10[i])

    kappa3=rep(0, 2)
    for(k in 1:n1)
    { if(y10[k]<=y10[i])
      kappa3=kappa3+apply(score1[y10>=y10[k],,drop=F], 2, sum)*(1-d10[k])/(sum(y10>=y10[k]))^2
    }

    kappa.arm1[i,]=kappa1+kappa2-kappa3
    }


    kappa.arm0=matrix(0, n0, 2)
    for(i in 1:n0)
    {kappa1=score0[i,]

    kappa2=apply(score0[y00>=y00[i],,drop=F], 2, sum)*(1-d00[i])/sum(y00>=y00[i])

    kappa3=rep(0, 2)
    for(k in 1:n0)
    {if(y00[k]<=y00[i])
      kappa3=kappa3+apply(score0[y00>=y00[k],,drop=F], 2, sum)*(1-d00[k])/(sum(y00>=y00[k]))^2
    }

    kappa.arm0[i,]=kappa1+kappa2-kappa3
    }


    if(type=="difference")
    {
      A=cbind(c(n1+n0, n1), c(n1, n1))
      betainf=solve(A)%*%t(rbind(kappa.arm1, kappa.arm0))
    }

    if(type=="ratio" || type=="lossratio")
    {mu1=exp(beta0[1]+beta0[2])
    mu0=exp(beta0[1])
    A=cbind(c(n1*mu1+n0*mu0, n1*mu1), c(n1*mu1, n1*mu1))

    betainf=solve(A)%*%t(rbind(kappa.arm1, kappa.arm0))
    }


    aug=rbind(x1*(1-pi), -x0*pi)

    fit=lm(t(betainf)~aug-1)

    if(type=="difference")
    {beta0=beta0
    se0=sqrt(diag(betainf%*%t(betainf)))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    cilow=beta0-se0*abs(qnorm((1-conf.int)/2))
    cihigh=beta0+se0*abs(qnorm((1-conf.int)/2))
    result.ini=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, label.low=cilow, label.upp=cihigh)
    colnames(result.ini)=c("coef", "se(coef)", "z", "p", label.low, label.upp)

    beta.aug=beta0-apply(aug%*%fit$coef,2,sum)
    se.aug=sqrt(diag(t(fit$res)%*%fit$res))*sqrt((n1+n0)/(n1+n0-p))
    z.aug=beta.aug/se.aug
    p.aug=1-pchisq(z.aug^2, 1)
    cilow.aug=beta.aug-se.aug*abs(qnorm((1-conf.int)/2))
    cihigh.aug=beta.aug+se.aug*abs(qnorm((1-conf.int)/2))
    result.aug=cbind(coef=beta.aug, "se(coef)"=se.aug, z=z.aug, p=p.aug, label.low=cilow.aug, label.upp=cihigh.aug)
    colnames(result.aug)=c("coef", "se(coef)", "z", "p", label.low, label.upp)
    }

    if(type=="ratio" ||  type=="lossratio")
    {beta0=beta0
    se0=sqrt(diag(betainf%*%t(betainf)))
    z0=beta0/se0
    p0=1-pchisq(z0^2, 1)
    cilow=beta0-se0*abs(qnorm((1-conf.int)/2))
    cihigh=beta0+se0*abs(qnorm((1-conf.int)/2))
    result.ini=cbind(coef=beta0, "se(coef)"=se0, z=z0, p=p0, "exp(coef)"=exp(beta0), label.low=exp(cilow), label.upp=exp(cihigh))
    colnames(result.ini)=c("coef", "se(coef)", "z", "p", "exp(coef)", label.low, label.upp)

    beta.aug=beta0-apply(aug%*%fit$coef,2,sum)
    se.aug=sqrt(diag(t(fit$res)%*%fit$res))*sqrt((n1+n0)/(n1+n0-p))
    z.aug=beta.aug/se.aug
    p.aug=1-pchisq(z.aug^2, 1)
    cilow.aug=beta.aug-se.aug*abs(qnorm((1-conf.int)/2))
    cihigh.aug=beta.aug+se.aug*abs(qnorm((1-conf.int)/2))
    result.aug=cbind(coef=beta.aug, "se(coef)"=se.aug, z=z.aug, p=p.aug, "exp(coef)"=exp(beta.aug), label.low=exp(cilow.aug), label.upp=exp(cihigh.aug))
    colnames(result.aug)=c("coef", "se(coef)", "z", "p", "exp(coef)", label.low, label.upp)

    }

    rownames(result.ini)=rownames(result.aug)=c("intercept", "arm")

    return(list(result.ini=result.ini, result.aug=result.aug))
  }
}


##################################################################################

# #
# data=read.table("c:/temp/data-e4a03-covs.csv", sep=",", head=T, na.string="NA")
# data=data[,-1]
# data$stage[is.na(data$stage)]=2
# y=data$time
# delta=data$status

# arm=data[,3]
# tau0=40

# x=as.matrix(data[,c(3, 4,7,8)])
# rmstreg(y, delta, x, arm, tau=tau0, type="difference")
# rmstreg(y, delta, x, arm, tau=tau0, type="ratio")
# rmstreg(y, delta, x, arm, tau=tau0, type="lossratio")

# x=as.matrix(data[, c(4,7,8)])
# rmstaug(y, delta, x, arm, tau=tau0, type="difference")
# rmstaug(y, delta, x, arm, tau=tau0, type="ratio")
# rmstaug(y, delta, x, arm, tau=tau0, type="lossratio")


# x=arm
# rmstreg(y, delta, x, arm, tau=tau0, type="difference")
# rmstreg(y, delta, x, arm, tau=tau0, type="ratio")
# rmstreg(y, delta, x, arm, tau=tau0, type="lossratio")

# x=data[,4]
# rmstaug(y, delta, x, arm, tau=tau0, type="difference")
# rmstaug(y, delta, x, arm, tau=tau0, type="ratio")
# rmstaug(y, delta, x, arm, tau=tau0, type="lossratio")
NULL


#' @name surv2sample
#' @aliases surv2sample
#' @title Inference of model-free between-group contrasts with censored survival data
#' @description Performs inference of several model-free group contrast measures, which include difference/ratio of cumulative incidence rates, quantiles, restricted mean survival times (RMST), and integrated survival rates.
#' @usage  surv2sample(time, status, arm, npert=1000,
#'                     timepoints=c(12, 24, 36, 40), quanprobs=c(0.1, 0.15, 0.2),
#'                     tau_start=0, tau, SEED=NULL, procedure="KM", conf.int=0.95)
#' @param time The follow-up time.
#' @param status The censoring indicator, 1=event, and 0=censoring.
#' @param arm The indicator for groups to compare 1/0.
#' @param npert The number of resampling. The default is 1000.
#' @param timepoints specifies the time points at which difference and ratio of the survival rates are computed.
#' @param quanprobs specifies the probabilities at which difference and ratio of the corresponding quantiles are computed.
#' @param tau_start The value indicates time point on the follow-up time to calculate the restricted mean survival time beyond the time point. The default is 0.
#' @param tau The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time. (i.e., the minimum of the largest observed time in each of the two groups)
#' @param SEED A random seed used for the resampling. Default is NULL.
#' @param procedure Specifies the inference procedure. A non-parametric procedure by the method of Kaplan-Meier ("KM") is the default. Another option is a parametric inference procedure by fitting a generalized gamma distribution to each group ("GG").
#' @param conf.int The level for computation of the confidence intervals. The default is 0.95.
#' @author Hajime Uno, Miki Horiguchi
#' @references
#' Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.
#'
#' Zhao L, Tian L, Uno H, Solomon S, Pfeffer M, Schindler J, Wei LJ. Utilizing the integrated difference of two survival functions to quantify the treatment contrast for designing, monitoring, and analyzing a comparative clinical study. Clinical Trials 2012, 9, 570-577.
#' @examples
#' D=pbc.sample()
#' surv2sample(D$time, D$status, D$group, npert=500, timepoints=c(2,4,6,8),
#' quanprobs =c(0.2, 0.3), tau=8, procedure="KM")
NULL


#####################################
# KM2.pert --hidden
# Perturbation: ver003 add averages
#####################################
# Without rectangle

KM2.pert <- function(time, status, npert=300, timepoints=c(12,24,36,40), quanprobs=c(0.1, 0.15, 0.2), tau_start=0, tau){

  indata=cbind(time, status)
  N=nrow(indata)
  KK=npert+1
  k.time=length(timepoints) ; p.time=matrix(0, nrow=KK, ncol=k.time)
  k.quan=length(quanprobs)  ; q.time=matrix(0, nrow=KK, ncol=k.quan)
  rmst=rep(0, KK)
  p.time.ave =rep(0, KK)
  q.time.ave =rep(0, KK)

  #if(is.null(tau)) tau=max(time)   ####changed

  #===========================================
  # observed
  #===========================================
  i=1
  ft= survfit(Surv(indata[,1], indata[,2])~1)

  #--- restricted mean survival time ---
  if(tau_start!=0){
    rmst[i]=summary(ft, rmean=tau)$table[5] - summary(ft, rmean=tau_start)$table[5]
  }else{
    rmst[i]=summary(ft, rmean=tau)$table[5]
  }


  #--- t-year survival ---
  for (k in 1:k.time){
    idx=ft$time<timepoints[k] ; p.time[i, k]=min(ft$surv[idx])}

  #--- quantiles ----------
  for (k in 1:k.quan){

    ######### ver004 #########

    idx1=round(ft$surv,digits=14) <= (1-quanprobs[k]) ;
    idx2=round(ft$surv,digits=14) <  (1-quanprobs[k]) ;
    if(sum(as.numeric(idx1))==0){
      q.time[i, k]=NA
    }else{
      if(sum(as.numeric(idx1)) == sum(as.numeric(idx2))){
        q.time[i, k]= min(ft$time[idx1])
      }else{
        q.time[i, k]=(min(ft$time[idx1]) + min(ft$time[idx2]))/2
      }
    }
  }
  #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

  #===========================================
  # perturbation
  #===========================================
  for (i in 2:KK){

    ft= survfit(Surv(indata[,1], indata[,2])~1, weight=rexp(N))

    #--- restricted mean survival time ---
    if(tau_start!=0){
      rmst[i]=summary(ft, rmean=tau)$table[5] - summary(ft, rmean=tau_start)$table[5]   ###changed
    }else{
      rmst[i]=summary(ft, rmean=tau)$table[5]
    }

    #--- t-year survival ---
    for (k in 1:k.time){
      idx=ft$time<timepoints[k] ; p.time[i, k]=min(ft$surv[idx])}

    #--- quantiles ----------
    for (k in 1:k.quan){

      ######### ver004 #########

      idx1=round(ft$surv,digits=14) <= (1-quanprobs[k]) ;
      idx2=round(ft$surv,digits=14) <  (1-quanprobs[k]) ;
      if(sum(as.numeric(idx1)) == sum(as.numeric(idx2))){
        q.time[i, k]= min(ft$time[idx1])
      }else{
        q.time[i, k]=(min(ft$time[idx1]) + min(ft$time[idx2]))/2
      }
    }

    #--- Average of t-year survivals and average percentiles (ver003) ---
    p.time.ave[i]=mean(p.time[i,])
    q.time.ave[i]=mean(q.time[i,])

  }
  #===========================================


  #--- output ---
  Z=list()
  Z$percentiles = data.frame(q.time); colnames(Z$percentiles)=quanprobs
  Z$tyearprobs  = data.frame(p.time) ; colnames(Z$tyearprobs) = timepoints
  Z$rmst        = rmst
  Z$tau         = tau
  Z$tau_start   = tau_start

  Z$tyearprobs.ave=p.time.ave
  Z$percentiles.ave=q.time.ave

  return(Z)
}
NULL


#############################################
# GG2.boot  --hidden
# GG2-02-boot-ver002.R: ver003 add averages
#############################################
GG2.boot=function(time, status, npert=300, timepoints=c(12,24,36,40), quanprobs=c(0.1, 0.15, 0.2), tau){

  indata=cbind(time, status)
  N=nrow(indata)
  KK=npert+1
  k.time=length(timepoints) ; p.time=matrix(0, nrow=KK, ncol=k.time)
  k.quan=length(quanprobs)  ; q.time=matrix(0, nrow=KK, ncol=k.quan)
  rmst=rep(0, KK)
  p.time.ave =rep(0, KK)
  q.time.ave =rep(0, KK)

  #if(is.null(tau)) tau=max(time[status==1])


  #===========================================
  # observed
  #===========================================
  i=1
  ft=flexsurvreg(Surv(indata[,1], indata[,2])~1, dist="gengamma")
  parm=ft$res[,1]

  #--- restricted mean survival time ---
  integrand=function(x){1-pgengamma(x, mu=parm[1], sigma = parm[2], Q=parm[3])}
  aa=integrate(integrand, lower=0, upper=tau)
  rmst[i]=aa$value

  #--- t-year survival ---
  for (k in 1:k.time){
    p.time[i, k]=1-pgengamma(timepoints[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

  #--- quantiles ----------
  for (k in 1:k.quan){
    q.time[i, k]=qgengamma(quanprobs[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

  #--- Average of t-year survivals and average percentiles (ver003) ---
  p.time.ave[i]=mean(p.time[i,])
  q.time.ave[i]=mean(q.time[i,])

  #===========================================
  # bootstrap
  #===========================================
  for (i in 2:KK){

    bs=sample(N, replace=TRUE)
    bdata=indata[bs,]

    ft=flexsurvreg(Surv(bdata[,1], bdata[,2])~1, dist="gengamma")
    parm=ft$res[,1]

    #--- restricted mean survival time ---
    integrand=function(x){1-pgengamma(x, mu=parm[1], sigma = parm[2], Q=parm[3])}
    aa=integrate(integrand, lower=0, upper=tau)
    rmst[i]=aa$value

    #--- t-year survival ---
    for (k in 1:k.time){
      p.time[i, k]=1-pgengamma(timepoints[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

    #--- quantiles ----------
    for (k in 1:k.quan){
      q.time[i, k]=qgengamma(quanprobs[k], mu=parm[1], sigma = parm[2], Q=parm[3])}

    #--- Average of t-year survivals and average percentiles (ver003) ---
    p.time.ave[i]=mean(p.time[i,])
    q.time.ave[i]=mean(q.time[i,])

  }
  #===========================================


  #--- output ---
  Z=list()
  Z$percentiles=data.frame(q.time); colnames(Z$percentiles)= quanprobs
  Z$tyearprobs=data.frame(p.time) ; colnames(Z$tyearprobs) = timepoints
  Z$rmst=rmst
  Z$tau=tau
  Z$tyearprobs.ave=p.time.ave
  Z$percentiles.ave=q.time.ave

  return(Z)
}
NULL







#' @export
###########################################################################
# surv2sampleComp (ver002) add average t-years and percentiles
# surv2sampleComp (ver003) add rmst beyond tau_start (without rectangle)
###########################################################################
surv2sample <- function(time, status, arm, npert=1000, timepoints=c(12,24,36,40), quanprobs =c(0.1, 0.15, 0.2), tau_start=0, tau, SEED=NULL, procedure="KM", conf.int=0.95){

  if(!is.null(SEED)){ set.seed(SEED) }

  #==================================
  #  initial check
  #==================================
  #--- tau ---
  idx=arm==0; tt=time[idx]; tau0max=max(tt)
  idx=arm==1; tt=time[idx]; tau1max=max(tt)
  tau_max = min(tau0max, tau1max)

  #---------------------
  if(!is.null(tau)){
    if(tau > tau_max){
      stop(paste("The truncation time, tau, needs to be shorter than or equal to the minimum of the largest observed time on each of the two groups: ", round(tau_max, digits=3)))
    }
  }
  #---------------------
  if(is.null(tau)){
    stop(paste("The truncation time, tau, was not specified. It needs to be specified by users."))
  }
  #---------------------

  if(procedure=="KM"){
    idx=arm==0; km0=KM2.pert(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau_start=tau_start, tau=tau)
    idx=arm==1; km1=KM2.pert(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau_start=tau_start, tau=tau)
  }

  if(procedure=="GG"){
    idx=arm==0; km0=GG2.boot(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
    idx=arm==1; km1=GG2.boot(time[idx], status[idx], npert=npert, timepoints=timepoints, quanprobs = quanprobs, tau=tau)
  }

  q_ci = function(x){quantile(x, prob=c((1-conf.int)/2, 1-(1-conf.int)/2))}

  time_interval = tau - tau_start #-- for calculating the integraated RMST --

  #--- each arm ---
  # wk0=cbind(km0$rmst, tau-km0$rmst, km0$tyearprobs, km0$percentiles)
  # wk1=cbind(km1$rmst, tau-km1$rmst, km1$tyearprobs, km1$percentiles)
  wk0=cbind(km0$rmst, tau-km0$rmst, km0$tyearprobs, km0$percentiles, km0$tyearprobs.ave, km0$percentiles.ave)
  wk1=cbind(km1$rmst, tau-km1$rmst, km1$tyearprobs, km1$percentiles, km1$tyearprobs.ave, km1$percentiles.ave)
  se0=apply(wk0[-1,], 2, sd)
  se1=apply(wk1[-1,], 2, sd)
  ci0=apply(wk0[-1,], 2, q_ci)
  ci1=apply(wk1[-1,], 2, q_ci)

  # measure=c("RMST","Loss time", paste("Prob at",timepoints), paste("Quantile at", quanprobs*100,"%"))
  measure=c("RMST","Loss time", paste("Prob at",timepoints), paste("Quantile at", quanprobs*100,"%"), "Ave of t-year event rates","Ave percentiles")

  out.group0=cbind(t(wk0[1,]), t(ci0), se0)
  out.group1=cbind(t(wk1[1,]), t(ci1), se1)
  rownames(out.group0)=measure
  rownames(out.group1)=measure
  colnames(out.group0)=c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"), "SE")
  colnames(out.group1)=c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"), "SE")

  #--- contrast ---
  K=ncol(wk0)
  outq=c()
  for (k in 1:K){
    q0=wk0[,k] ; q1=wk1[,k]
    diff01=q0-q1
    diff10=q1-q0
    ratio01=q0/q1
    ratio10=q1/q0
    outq=cbind(outq, diff01, diff10, ratio01, ratio10)
  }

  measures=rep(measure, each=4)
  contrast=rep(c("Group0-Group1","Group1-Group0","Group0/Group1","Group1/Group0"), K)

  #======================================
  if(procedure=="KM"){


    #--- p-val --
    se=apply(outq[-1,], 2, sd)
    pval=pnorm(-abs(outq[1,])/se)*2

    idx=contrast=="Group0/Group1"|contrast=="Group1/Group0"
    for (j in 1:length(pval)){
      if(contrast[j]=="Group0/Group1" | contrast[j]=="Group1/Group0"){
        for (m in 1:length(measure)){
          if(measures[j]==measure[m]){
            se[j]=sqrt( (se0[m]/out.group0[m,1])^2 + (se1[m]/out.group1[m,1])^2)
          }
        }
      }
    }
    pval[idx]=pnorm(-abs(log(outq[1,idx]))/se[idx])*2

    #--- confidence interals --
    lower=outq[1,] - abs(qnorm((1-conf.int)/2))*se
    upper=outq[1,] + abs(qnorm((1-conf.int)/2))*se
    idx=contrast=="Group0/Group1"|contrast=="Group1/Group0"
    lower[idx]=exp(log(outq[1,idx]) - abs(qnorm((1-conf.int)/2))*se[idx])
    upper[idx]=exp(log(outq[1,idx]) + abs(qnorm((1-conf.int)/2))*se[idx])



    #--- results ---
    out.contrast=cbind(outq[1,], lower, upper, pval)
    rownames(out.contrast)=paste(measures, contrast)
    colnames(out.contrast)=c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"), "p-val")
    inf.method="Perturbation resampling"

  }
  #======================================
  if(procedure=="GG"){

    #--- confidence interals (bootstrap percentile) --
    cband=apply(outq[-1,], 2, q_ci)
    lower=cband[1,]
    upper=cband[2,]

    out.contrast=cbind(outq[1,], lower, upper)
    rownames(out.contrast)=paste(measures, contrast)
    colnames(out.contrast)=c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"))
    inf.method="Bootstrap percentile method"

  }
  #======================================


  #--- intergrated rmst ---
  int.surv.group0 = cbind(out.group0[1,1]/time_interval, out.group0[1,2]/time_interval, out.group0[1,3]/time_interval)
  int.surv.group1 = cbind(out.group1[1,1]/time_interval, out.group1[1,2]/time_interval, out.group1[1,3]/time_interval)

  int.surv.diff01 = cbind(out.contrast[contrast=="Group0-Group1",][1,1]/time_interval, out.contrast[contrast=="Group0-Group1",][1,2]/time_interval, out.contrast[contrast=="Group0-Group1",][1,3]/time_interval, out.contrast[contrast=="Group0-Group1",][1,4])
  int.surv.diff10 = cbind(out.contrast[contrast=="Group1-Group0",][1,1]/time_interval, out.contrast[contrast=="Group1-Group0",][1,2]/time_interval, out.contrast[contrast=="Group1-Group0",][1,3]/time_interval, out.contrast[contrast=="Group1-Group0",][1,4])

  out.int.surv        = rbind(int.surv.group0, int.surv.group1)
  out.int.surv.diff   = rbind(int.surv.diff01, int.surv.diff10)

  rownames(out.int.surv) = c("Integrated survival rate Group0", "Integrated survival rate Group1")
  colnames(out.int.surv) = c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"))

  rownames(out.int.surv.diff) = c("Integrated survival rate Group0-Group1", "Integrated survival rate Group1-Group0")
  colnames(out.int.surv.diff) = c("Est.", paste0("Lower ", conf.int*100, "%"), paste0("Upper ", conf.int*100, "%"), "p-val")


  #--- output ---
  Z=list()
  Z$procedure            = procedure
  Z$method               = inf.method
  Z$survfit              = survfit(Surv(time, status)~arm)
  Z$tau_start            = tau_start
  Z$tau                  = tau
  Z$time_interval        = time_interval
  Z$npert                = npert
  Z$timepoints           = timepoints
  Z$quanprobs            = quanprobs
  Z$contrast.all         = out.contrast
  Z$group0               = out.group0
  Z$group1               = out.group1
  Z$RMST                 = out.contrast[measures=="RMST",]
  Z$RMLT                 = out.contrast[measures=="Loss time",]
  Z$contrast.diff10      = out.contrast[contrast=="Group1-Group0",]
  Z$contrast.diff01      = out.contrast[contrast=="Group0-Group1",]
  Z$contrast.ratio01     = out.contrast[contrast=="Group0/Group1",]
  Z$contrast.ratio10     = out.contrast[contrast=="Group1/Group0",]
  Z$integrated_surv      = out.int.surv
  Z$integrated_surv.diff = out.int.surv.diff

  class(Z)="surv2sample"
  return(Z)

}
NULL
