\name{plot.surv2sample}
\alias{plot.surv2sample}
\title{Plot method for surv2sample objects}

\description{
Creates plots from a surv2sample object. 
}

\usage{
\method{plot}{surv2sample}(x, measure=NULL, baseline=0, \dots)
}

\arguments{
  \item{x}{
\code{surv2sample} object
}
  \item{measure}{
The type of measure used for the plot. When default(=NULL), \code{plot.survfit()} is called and KM plots are given.
When "relative time" is specified, a plot of relative percentiles with corresponidng 0.95 confidence intervals is generatead.  
}
  \item{baseline}{
Indicates the baseline group, 0/1. Default is 0.
}
  \item{\dots}{
for future method
}
}

\seealso{\code{plotCI}(plotrix)}


