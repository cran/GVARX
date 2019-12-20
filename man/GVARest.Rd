\name{GVARest}
\alias{GVARest}

\title{Estimate country-specific VAR in a GVAR setting
}

\encoding{latin1}

\description{Estimate country-specific VAR in a GVAR setting}

\usage{GVARest(data,p,FLag, lag.max, type="const", ic,weight.matrix=NULL)
}

\arguments{
  \item{data}{Dataframe is a strictly balanced panel data format,the first column is cross-section ID,and the second column is Time. For the sake of identification, both columns must be named by, respectively, id and Time.
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

\value{

  \item{gvar}{Country-specific GVAR output list}
  \item{White}{Coefficient estimates with White robust covariance}
  \item{NWHAC}{Coefficient estimates withNewy-West robust covariance}
  \item{p}{Number of lags for endogeneous variables in VAR}
  \item{K}{Number of lags for Ft variables in VAR}
  \item{type}{Model specificaiton. As in package vars, we have four selection: "none","const","trend", and "both".}
  \item{datamat}{input data=data}
  \item{lagmatrix}{GVAR's Country-secific optimal lag number.}
  \item{lagmatrix1}{VAR's Country-secific optimal lag number.}
  \item{exoLag}{Ft lags}
  \item{Ft}{Foreign variables}
  \item{NAMES}{Names of countries}
  \item{gvarRSD}{Country-specific GVAR residuals}
  \item{varRSD}{VAR residuals}
  \item{weight}{weight.matrix}
}

\references{Mauro Filippo di and Pesaran H. M. (2013) The GVAR Handbook-- Structure and Applications of a Macro Model of the Global Economy for Policy. Oxford University Press.
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
mainOUTPUT = GVARest(data=PriceVol,p,FLag,lag.max,type,ic,weight.matrix)

mainOUTPUT$lagmatrix    # Country-specific GVAR lags
mainOUTPUT$gvar
mainOUTPUT$gvar[[1]]
coef(mainOUTPUT$gvar[[17]])
mainOUTPUT$White[[17]]
mainOUTPUT$NWHAC[[17]][1]

}
