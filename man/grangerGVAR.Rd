\name{grangerGVAR}
\alias{grangerGVAR}

\title{
Granger Causality for bivariate vector autoregression model augmented by foreign variables
}

\encoding{latin1}

\description{
Granger Causality for bivariate vector autoregression (VAR) model augmented by foreign variables
}
\usage{
grangerGVAR(data,p=2,FLag=2, type="const", lag.max,  ic,weight.matrix=NULL)
}

\arguments{
  \item{data}{Dataframe is a strictly balanced panel data format, the first column is cross-section ID, and the second column is timeID. For the sake of identification, both columns must be named by, respectively, ID and timeID.
}

 \item{p}{The number of lag for Xt matrix.
}

 \item{FLag}{The number of lag for foreign variables in country-specific VAR

}

 \item{lag.max}{The maximal number of lag for estimating country-specific VAR
}

 \item{type}{Model specificaiton for VAR. As in package vars, we have four selection: "none","const","trend", "both".
}

 \item{ic}{Information criteria for optimal lag.As in package vars, we have four selection: "AIC", "HQ", "SC", and "FPE".
}

 \item{weight.matrix}{Bilateral trade weight matrix for computing foreign variables. If the computation of foreign variables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign variables are weighted on a year-to-year basis, then weight.matrix must be a "list, with the same length as the weighting frequency.
}
}

\details{
This functions implements country-specific Granger causality test for bivariate VAR only. We restrict it to the case of two-variable and country specific, becasue Granger causality will be a little complex when the VAR system has more than two variables. In a global VAR case, the reduced large matrix Xt may not be suitable for Granger causality.
}

\value{
  \item{y1GCy2.var}{VAR Granger causality: y1(The 1st variable) does not Granger cause y2 (The 2nd variable)}
  \item{y2GCy1.var}{VAR Granger causality: y2(The 2nd variable) does not Granger cause y1 (The 1st variable)}
  \item{y1GCy2.gvar}{GVAR Granger causality: y1(The 1st variable) does not Granger cause y2 (The 2nd variable)}
  \item{y2GCy1.gvar}{GVAR Granger causality: y2(The 2nd variable) does not Granger cause y1 (The 1st variable)}


}
\references{Mauro Filippo di and Pesaran H. M. (2013) The GVAR Handbook-- Structure and Applications of a Macro Model of the Global Economy for Policy. Oxford University Press.

}
\author{
Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\examples{
data("PriceVol")
data("tradeweight1")
data("tradeweightx")

p=2
FLag=2
type="const"
lag.max=15
ic="SC"
weight.matrix=tradeweightx

GC_OUTPUT = grangerGVAR(data=PriceVol, p, FLag, type, lag.max, ic, weight.matrix)

# VAR:y1(Price) does NOT Granger Casuse y2 (logVol)
round(GC_OUTPUT$y1GCy2.var, 4)

# VAR: y2(logVol) does NOT Granger Casuse y1 (Price)
round(GC_OUTPUT$y2GCy1.var, 4)

# GVAR: y1(Price) does NOT Granger Casuse y2 (logVol)
round(GC_OUTPUT$y1GCy2.gvar, 4)

# GVAR: y2(logVol) does NOT Granger Casuse y1 (Price)
round(GC_OUTPUT$y2GCy1.gvar, 4)


}

