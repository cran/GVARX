\name{GVECMest}
\alias{GVECMest}

\title{Estimate country-specific Engle-Granger VECM in a Global VECM setting
}

\encoding{latin1}

\description{Estimate country-specific Engle-Granger VECM in a Global VECM setting}

\usage{GVECMest(data,p=2,lag.max=NULL, type="const", ic,weight.matrix=NULL)
}

\arguments{
  \item{data}{
  Dataframe is a strictly balanced panel data format,the first column is cross-section ID,and the second column is Time. For the sake of identification, both columns must be named by, respectively, id and Time.
}

 \item{p}{The number of lag for Xt matrix. Foreign variables are set by FLag=p+1. Current version restricts p <= 2 for simplicity, which aims at avoiding too many paramaters in low-frequency data of many variables and many countries. It will be relaxed soon.
}


 \item{lag.max}{The maximal number of lag for estimating country-specific VECM
}

 \item{type}{Model specificaiton for VECM. As in package VECMs, we have four selection: "none","const","trend", "both".
}

 \item{ic}{Information criteria for optimal lag.As in package VECMs, we have four selection: "AIC", "HQ", "SC", and "FPE".
}

 \item{weight.matrix}{Bilateral trade weight matrix for computing foreign VECMiables. If the computation of foreign VECMiables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign VECMiables are weighted on a year-to-year basis, then weight.matrix must be a "list, with the same length as the weighting frequency.
}


}

\value{

  \item{gvecm}{Country-specific GVECM output list}
  \item{White}{Coefficient estimates with White robust coVECMiance}
  \item{NWHAC}{Coefficient estimates withNewy-West robust coVECMiance}
  \item{p}{Number of lags for endogeneous VECMiables in VECM}
  \item{K}{Number of lags for Ft VECMiables in VECM}
  \item{type}{Model specificaiton. As in package VECMs, we have four selection: "none","const","trend", and "both".}
  \item{datamat}{input data=data}
  \item{lagmatrix}{GVECM's Country-secific optimal lag number.}
  \item{lagmatrix1}{VECM's Country-secific optimal lag number.}
  \item{exoLag}{Ft lags}
  \item{Ft}{Foreign VECMiables}
  \item{NAMES}{Names of countries}
  \item{gvecmRSD}{Country-specific Global VECM residuals}
  \item{vecmRSD}{VECM residuals}

}

\references{Mauro Filippo di and Pesaran H. M. (2013) The GVECM Handbook-- Structure and Applications of a Macro Model of the Global Economy for Policy. Oxford University Press.
}

\author{Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\examples{
data("PriceVol")
data("tradeweight1")
data("tradeweightx")

p=2
FLag=2
lag.max=15
type="const"
ic="SC"
weight.matrix=tradeweightx
mainOUTPUT = GVECMest(data=PriceVol,p,lag.max,type,ic,weight.matrix)

mainOUTPUT$lagmatrix    # Country-specific GVECM lags
mainOUTPUT$gvecm
mainOUTPUT$gvecm[[1]]
coef(mainOUTPUT$gvecm[[17]])
mainOUTPUT$White[[17]]
mainOUTPUT$NWHAC[[17]][1]

}
