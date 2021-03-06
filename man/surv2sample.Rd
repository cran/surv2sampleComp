% Generated by roxygen2: do not edit by hand
% Please edit documentation in R/surv2sample-ver003.R
\name{surv2sample}
\alias{surv2sample}
\title{Inference of model-free between-group contrasts with censored survival data}
\usage{
surv2sample(time, status, arm, npert=1000,
                    timepoints=c(12, 24, 36, 40), quanprobs=c(0.1, 0.15, 0.2),
                    tau_start=0, tau, SEED=NULL, procedure="KM", conf.int=0.95)
}
\arguments{
\item{time}{The follow-up time.}

\item{status}{The censoring indicator, 1=event, and 0=censoring.}

\item{arm}{The indicator for groups to compare 1/0.}

\item{npert}{The number of resampling. The default is 1000.}

\item{timepoints}{specifies the time points at which difference and ratio of the survival rates are computed.}

\item{quanprobs}{specifies the probabilities at which difference and ratio of the corresponding quantiles are computed.}

\item{tau_start}{The value indicates time point on the follow-up time to calculate the restricted mean survival time beyond the time point. The default is 0.}

\item{tau}{The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time. (i.e., the minimum of the largest observed time in each of the two groups)}

\item{SEED}{A random seed used for the resampling. Default is NULL.}

\item{procedure}{Specifies the inference procedure. A non-parametric procedure by the method of Kaplan-Meier ("KM") is the default. Another option is a parametric inference procedure by fitting a generalized gamma distribution to each group ("GG").}

\item{conf.int}{The level for computation of the confidence intervals. The default is 0.95.}
}
\description{
Performs inference of several model-free group contrast measures, which include difference/ratio of cumulative incidence rates, quantiles, restricted mean survival times (RMST), and integrated survival rates.
}
\examples{
D=pbc.sample()
surv2sample(D$time, D$status, D$group, npert=500, timepoints=c(2,4,6,8),
quanprobs =c(0.2, 0.3), tau=8, procedure="KM")
}
\references{
Tian L, Zhao L, Wei LJ. Predicting the restricted mean event time with the subject's baseline covariates in survival analysis. Biostatistics 2014, 15, 222-233.

Zhao L, Tian L, Uno H, Solomon S, Pfeffer M, Schindler J, Wei LJ. Utilizing the integrated difference of two survival functions to quantify the treatment contrast for designing, monitoring, and analyzing a comparative clinical study. Clinical Trials 2012, 9, 570-577.
}
\author{
Hajime Uno, Miki Horiguchi
}
