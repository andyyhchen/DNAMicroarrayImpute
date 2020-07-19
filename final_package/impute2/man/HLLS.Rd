\name{HLLS}
\alias{HLLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
HLLS DNA imputation
}
\description{
	Hybird Local Least Square imputation of DNA microarry missing value 
}
\usage{
HLLS(xmiss, K = 15, sim.method = "EuDist", Niter = 2, order = 1)
}
\arguments{
  \item{xmiss}{
An matrix with missing value
}
  \item{K}{
number of neighorhood selected to estimate missing value
}
  \item{sim.method}{
method to select neighborhood
}
  \item{Niter}{
times of iteration
}

  \item{order}{
just order!
}
}
\author{
}
\note{
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##---- Should be DIRECTLY executable !! ----
##-- ==>  Define data, use random,
##--	or do  help(data=index)  for the standard data sets.

## The function is currently defined as
function (xmiss, K = 15, sim.method = "EuDist", Niter = 2, order = 1) 
{
    xcomplete <- WSLLS(xmiss, K, sim.method, order)
    miss.gene <- is.na(xmiss)
    miss.row <- which(rowSums(miss.gene) != 0)
    for (h in 1:Niter) {
        xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, 
            function(row) {
                row.miss <- which(is.na(row))
                sim <- similarityCal(row[-row.miss], xcomplete[, 
                  -row.miss], sim.method)
                sim.id <- order(sim, decreasing = T)[2:K + 1]
                x.diag <- diag(sim[sim.id]^order)
                row[row.miss] <- t(x.diag \%*\% xcomplete[sim.id, 
                  row.miss, drop = FALSE]) \%*\% ginv(t(x.diag \%*\% 
                  xcomplete[sim.id, -row.miss, drop = FALSE])) \%*\% 
                  row[-row.miss, drop = FALSE]
                return(row)
            }))
    }
    return(xcomplete)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
