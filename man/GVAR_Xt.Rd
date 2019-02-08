\name{GVAR_Xt}
\alias{GVAR_Xt}

\title{Compute the G0, G1, G2, and F1, F2 matrices for filtering Xt
}

\encoding{latin1}

\description{Compute the G0, G1, G2, and F1, F2 matrices for filtering Xt}

\usage{GVAR_Xt(data,p=1,FLag, lag.max, type="const", ic, weight.matrix=NULL)
}

\arguments{
  \item{data}{
Dataframe is a strictly balanced panel data format, the first column is cross-section ID, and the second column is timeID. For the sake of identification, both columns must be named by, respectively, ID and timeID.
}

 \item{p}{The number of lag for Xt matrix.
}

 \item{FLag}{The number of lag for foreign variables in country-specific VAR

}

 \item{lag.max}{The maximal number of lag for estimating country-specific VAR
}

 \item{type}{Model specificaiton for VAR. As in package vars, we have four selection: "none","const","trend", "both".
}

 \item{ic}{
Information criteria for optimal lag.As in package vars, we have four selection: "AIC", "HQ", "SC", "FPE".
}

 \item{weight.matrix}{Bilateral trade weight matrix for computing foreign variables. If the computation of foreign variables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign variables are weighted on a year-to-year basis, then weight.matrix must be a "list", with the same length as the weighting frequency.
}


}

\details{This function generates several matrices of Eq.(2.6) in Filippo and Pesaran(2013, P.17), which is resuired to recursively filter Xt; besides, it also re-calculates the transformed residuals. In this version, we do not include the impulse responses function(IRF), because the IRF can be computed by these matrices and residuals easily. We will not update it until the next version.
}

\value{

  \item{G0}{Matrix G0 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{G1}{Matrix G1 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{G2}{Matrix G2 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{F1}{Matrix F1 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{F2}{Matrix F2 of Eq.(2.6) in Filippo and Pesaran(2013, P.17)}
  \item{lagmatrix}{Country-secific optimal lag number.}
  \item{newRESID}{New residuals=epsilon in Filippo and Pesaran (2013, P.17)}
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
FLag=2
lag.max=15
type="const"
ic="SC"
weight.matrix=tradeweightx

Result=GVAR_Xt(data=PriceVol,p,FLag,lag.max,type,ic, weight.matrix)
Result$G0
Result$G1
#Result$F1
Result$lagmatrix
Result$newRESID

}
