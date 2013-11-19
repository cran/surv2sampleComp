\name{rmstaug}
\alias{rmstaug}
\title{Adjusted difference/ratio of restricted mean survival times}
\description{
Compares restricted mean survival time between two groups, adjusting for imbalance of baseline factors. 
}
\usage{
rmstaug(y, delta, x, arm, tau, type = "difference")
}
\arguments{
  \item{y}{
The follow-up time.
}
  \item{delta}{
The censoring indicator, 1=event, and 0= censoring.
}
  \item{x}{
The covariate matrix. The group indicator, \code{arm} (below) should not be included in this matrix. 
}
  \item{arm}{
The group indicator, 1/0.
}
  \item{tau}{
The value indicates the restricted time point on the follow-up time to calculate the restricted mean survival time. 
}
  \item{type}{
The type of the between-group contrast measure: "difference"(default), "ratio" or "lossratio" 
} 
}
\references{
Tian, Lu; Zhao, Lihui; and Wei, LJ. "On the Restricted Mean Event Time in Survival Analysis." (February 2013). Harvard University Biostatistics Working Paper Series. Working Paper 156. http://biostats.bepress.com/harvardbiostat/paper156}
\author{
Lu Tian}
\examples{
D=pbc.sample()

rmstaug(D$time, D$status, D$covariates, D$group, tau=8, type="difference")

}
