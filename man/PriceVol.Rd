\name{PriceVol}
\alias{PriceVol}
\docType{data}

\title{Dataset price-volumn of 17 mareket indices
}

\encoding{latin1}

\description{A nine-year balanced panel price-volumn data of 17 mareket indices, 2006/8/30-2014/11/19
}
\usage{data("PriceVol")}
\format{
  A data frame with 0 observations on the following 2 variables.
  \describe{
    \item{\code{ID}}{Names of country, cross-section ID}
    \item{\code{Time}}{Time index}
    \item{\code{Ret}}{Daily returns computed by close-to-close}
    \item{\code{Vol}}{Daily transaction volumn, by log}
  }
}

\source{Yahoo finance
}

\examples{data(PriceVol)

}
