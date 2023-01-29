\name{GVAR_GF}
\alias{GVAR_GF}

\title{Compute the structural coefficients matrices G0, G1, G2, and F1, F2
}

\encoding{latin1}

\description{Compute the structural coefficients matrices G0, G1, G2, and F1, F2}

\usage{GVAR_GF(data,p, type="const",ic="AIC",weight.matrix)
}

\arguments{
  \item{data}{
Dataframe is a strictly balanced panel data format,the first column is cross-section ID,and the second column is Time. For the sake of identification, both columns must be named by, respectively, id and Time.
}

 \item{p}{The number of lag for Xt matrix. The number of lag for foreign variables in country-specific VAR FLag is set to be p+1. Current version restricts p <= 2 for simplicity, which aims at avoiding too many paramaters in low-frequency data of many variables and many countries. It will be relaxed soon.
}

 \item{type}{Model specificaiton for VAR. As in package vars, we have four selection: "none","const","trend", "both".
}

 \item{ic}{
Information criteria for optimal lag.As in package vars, we have four selection: "AIC", "HQ", "SC", "FPE".
}

 \item{weight.matrix}{Bilateral trade weight matrix for computing foreign variables. If the computation of foreign variables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign variables are weighted on a year-to-year basis, then weight.matrix must be a "list", with the same length as the weighting frequency.
}


}

\details{This function generates several structural coefficient matrices of Eq.(2.6) in Filippo and Pesaran(2013, P.17), which are required to compute IRF and multistep forecasts. Besides, it also re-calculates the transformed residuals. In this version, we do not include the impulse responses function(IRF), because the IRF can be computed by these matrices and residuals easily. We will not update it until the next version.
}

\value{

  \item{G0}{Matrix G0 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{G1}{Matrix G1 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{G2}{Matrix G2 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{F1}{Matrix F1 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{F2}{Matrix F2 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{lagmatrix}{Country-secific optimal lag number, which must be the same.}
  \item{RESID}{original residuals=u in Filippo and Pesaran (2013, P.17)}
  \item{newRESID}{New residuals=epsilon in Filippo and Pesaran (2013, P.17)}
  \item{fitted}{In-sample fitted values, or conditional mean}
  \item{data}{data used}
}


\references{Mauro Filippo di and Pesaran H. M. (2013) The GVAR Handbook-- Structure and Applications of a Macro Model of the Global Economy for Policy. Oxford University Press.
}

\author{Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}


\examples{
data("PriceVol")
data("tradeweightx")
data("tradeweight1")
p=2
type="const"
ic="SC"

Result=GVAR_GF(data=PriceVol,p,type,ic, weight.matrix=tradeweight1)
Result$G0
Result$G1
Result$G2
Result$F1
Result$F2
Result$lagmatrix
Result$RESID
Result$newRESID
Result$fitted
Result$data
#May use forecast::accuracy(Result$fitted[,1], Result$data[,1]) for performance.

}
