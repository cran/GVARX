\name{tradeweightx}
\alias{tradeweightx}
\docType{data}

\title{A nine-year bilateral trade weight matrix, 2006-2014
}

\encoding{latin1}


\description{A nine-year bilateral trade weight matrix, 2006-2014
}
\usage{data("tradeweightx")}
\format{
  A list with 17 by 17 matrix on the following variable.
  \describe{
    \item{\code{Australia}}{Bilateral trade weight matrix of Australia, 2014}
    \item{\code{Austria}}{Bilateral trade weight matrix of Austria, 2014}
    \item{\code{Belgium}}{Bilateral trade weight matrix of Belgium, 2014}
    \item{\code{Brazil}}{Bilateral trade weight matrix of Brazil, 2014}
    \item{\code{France}}{Bilateral trade weight matrix of France, 2014}
    \item{\code{UK}}{Bilateral trade weight matrix of UK, 2014}
    \item{\code{US}}{Bilateral trade weight matrix of US, 2014}
    \item{\code{Canada}}{Bilateral trade weight matrix of Canada, 2014}
    \item{\code{HongKong}}{Bilateral trade weight matrix of Hong Kong, 2014}
    \item{\code{Indonesia}}{Bilateral trade weight matrix of Indonesia, 2014}
    \item{\code{Malaysia}}{Bilateral trade weight matrix of Malaysia, 2014}
    \item{\code{Korea}}{Bilateral trade weight matrix of Korea, 2014}
    \item{\code{Mexico}}{Bilateral trade weight matrix of Mexico, 2014}
    \item{\code{Japan}}{Bilateral trade weight matrix of Japan, 2014}
    \item{\code{Swiss}}{Bilateral trade weight matrix of Swiss, 2014}
    \item{\code{China}}{Bilateral trade weight matrix of China, 2014}
    \item{\code{Taiwan}}{Bilateral trade weight matrix of Taiwan, 2014}

  }
}

\details{This example data is annual trade weight matrix, it is a list with length 9 (2006-2014).Each list is a year specific 17 by 17 trade weight matrix, the column names are 17 countries. Given column j, the row-wise elements are bilateral trade weights of country j. Make sure that the length of list must exactly match with the number of years. Because once you use this as tradewieght input matrix, R function will automatically compute foreign variables weighted year-by-year. Please make sure that the order of countries exactly matches the dataset's ID column.
}

\examples{
data(tradeweightx)
is.data.frame(tradeweightx)
}
