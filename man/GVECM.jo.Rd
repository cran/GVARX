\name{GVECM.jo}
\alias{GVECM.jo}

\title{Estimate country-specific Johansen test results in a Global VECM setting
}

\encoding{latin1}

\description{Estimate country-specific Johansen test results in a Global VECM setting}

\usage{
GVECM.jo(data,p=2,FLag=3,ecdet = "const", type = "eigen",spec = "longrun",
season = NULL,weight.matrix)
}

\arguments{
  \item{data}{
  Dataframe is a strictly balanced panel data format,the first column is cross-section ID,and the second column is Time. For the sake of identification, both columns must be named by, respectively, id and Time.
}

 \item{p}{The number of lag for Xt matrix.
}

 \item{FLag}{The number of lag for foreign VECMiables in country-specific VECM

}

 \item{ecdet}{Character, 'none' for no intercept in cointegration, 'const' for constant term in
cointegration and 'trend' for trend variable in cointegration.
}

 \item{type}{Model specificaiton for VECM. As in package VECMs, we have four selection: "none","const","trend", "both".
}

 \item{spec}{Determines the specification of the VECM, see details in pakcage urca.
}


 \item{season}{If seasonal dummies should be included, the data frequency must be set accordingly,i.e '4' for quarterly data.
}

 \item{weight.matrix}{Bilateral trade weight matrix for computing foreign VECMiables. If the computation of foreign VECMiables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign VECMiables are weighted on a year-to-year basis, then weight.matrix must be a "list, with the same length as the weighting frequency.
}


}

\value{

  \item{JO.test}{List object of country-specific Johansen test results}
  \item{VECMoutputs}{List object of country-specific VECM  results}
  \item{RESID}{List object of country-specific VECM  residuals, obtained by using vars::vec2var}

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
type="const"
ic="SC"
weight.matrix=tradeweight1
mainOUT.JO=GVECM.jo(data=PriceVol,p=2,FLag=3, weight.matrix=weight.matrix)
mainOUT.JO$JO.test

}
