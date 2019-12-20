\name{GVAR_Ft}
\alias{GVAR_Ft}

\title{Function to generate foreign variables
}

\encoding{latin1}

\description{Function to generate foreign variables
}

\usage{GVAR_Ft(data, weight.matrix=NULL)

}

\arguments{
  \item{data}{
  Dataframe is a strictly balanced panel data format,the first column is cross-section ID,and the second column is Time. For the sake of identification, both columns must be named by, respectively, id and Time.
}
 \item{weight.matrix}{
Bilateral trade weight matrix for computing foreign variables. If the computation of foreign variables are weighted by one weighting matrix, weight.matrix must be a "data.frame". If the computation of foreign variables are weighted on a year-to-year basis, then weight.matrix must be a "list", with the same length as the weighting frequency. If NULL, then it computes the foreign vriables by average.
}

}

\value{

 \item{Ft}{Weighted foerign variables as described in GVAR}

}


\references{Mauro Filippo di and Pesaran H. M. (2013) The GVAR Handbook-- Structure and Applications of a Macro Model of the Global Economy for Policy. Oxford University Press.
}

\author{Ho Tsung-wu <tsungwu@ntnu.edu.tw>, College of Management, National Taiwan Normal University.
}

\examples{
#=== Loading Data ===#
data("PriceVol")
data("tradeweight1")
data("tradeweightx")

#Generate country-specific foreign variables
Ft=GVAR_Ft(data=PriceVol,weight.matrix=tradeweight1)
k=17
head(Ft[[k]])
tail(Ft[[k]])
}
