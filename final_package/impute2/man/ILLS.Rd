\name{ILLS}
\alias{ILLS}
\title{
	Iterated Local Least Square impute
}
\description{
ILLS DNA Impute
}

\usage{
ILLS(xmiss, K = 15, sim.method = "EuDist", Niter = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmiss}{
	  matrix with missing value
}
  \item{K}{
	  number of neighbors
}
  \item{sim.method}{
	  method to select neighbors
}
  \item{Niter}{
	  times of iteration
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
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
function (xmiss, K = 15, sim.method = "EuDist", Niter = 2) 
{
    xcomplete <- RowAverage(xmiss)
    miss.gene <- is.na(xmiss)
    miss.row <- which(rowSums(miss.gene) != 0)
    for (h in 1:Niter) {
        xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, 
            function(row) {
                row.miss <- which(is.na(row))
                sim <- similarityCal(row[-row.miss], xcomplete[, 
                  -row.miss, drop = F], sim.method)
                sim.id <- order(sim, decreasing = T)[2:K + 1]
                row[row.miss] <- t(xcomplete[sim.id, row.miss, 
                  drop = FALSE]) \%*\% ginv(t(xcomplete[sim.id, 
                  -row.miss, drop = FALSE])) \%*\% row[-row.miss, 
                  drop = FALSE]
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
