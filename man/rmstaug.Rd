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
D=CompCase(pbc[,c(2:4,10:14)]) 
year=D$time/365.25
status=as.numeric(D[,2]==2)
arm=as.numeric(D$trt==2)
covs=as.matrix(D[,4:8])

rmstaug(year, status, covs, arm, tau=8, type="difference")

}
