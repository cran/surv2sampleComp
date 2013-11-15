\name{surv2sample}
\alias{surv2sample}
\title{Inference of model-free between-group contrasts with censored survival data
}
\description{
Performs inference of several model-free group contrast measures, which include difference/ratio of cumulative incidence rates, quantiles, and restricted mean survival times. 
}
\usage{
surv2sample(time, status, arm, npert = 1000, timepoints = c(12, 24, 36, 40), 
quanprobs = c(0.1, 0.15, 0.2), tau = NULL, SEED = NULL, procedure="KM")
}
\arguments{
  \item{time}{
The follow-up time.
}
  \item{status}{
The censoring indicator, 1=event, and 0= censoring.
}
  \item{arm}{
The indicator for groups to compare 1/0.
}
  \item{npert}{
The number of resampling. The default is 1000. 
}
  \item{timepoints}{
specifies the time points at which difference and ratio of the survival rates are computed.
}
  \item{quanprobs}{
specifies the probabilities at which difference and ratio of the corresponding quantiles are computed.
}
  \item{tau}{
The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time. If \code{NULL}, the function automatically calculates the largest time among non-censored observations for each group and picks the smaller one. 
}
  \item{SEED}{
A random seed used for the resampling. Default is \code{NULL}. 
}
  \item{procedure}{
Specifies the inference procedure. A non-parametric procedure by the method of Kaplan-Meier ("KM") is the default. Another option is a parametric inference procedure by fitting a generalized gamma distribution to each group ("GG").
}
}
\author{Hajime Uno}
\examples{
D=CompCase(pbc[,c(2:4,10:14)]) 
year=D$time/365.25
status=as.numeric(D[,2]==2)
arm=as.numeric(D$trt==2)

surv2sample(year, status, arm, npert=500, timepoints=c(2,4,6,8), 
                  quanprobs =c(0.25, 0.5), tau=8, procedure="KM")

}

