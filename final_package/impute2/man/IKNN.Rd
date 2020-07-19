\name{IKNN}
\alias{IKNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	IKNN DNA Imputation
}
\description{
Iterated K nearest neighbors imputation
}
\usage{
IKNN(xmiss, K = 15, sim.method = "EuDist", Niter = 2)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmiss}{
	  matrix with missing value
}
  \item{K}{
	  number of neighbors choosen
}

  \item{sim.method}{
	  method to select neigbors
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
    miss.exp <- lapply(miss.row, function(i) which(is.na(xmiss[i, 
        ])))
    for (h in 1:Niter) {
        xcomplete[miss.row, ] <- t(apply(xmiss[miss.row, ], 1, 
            function(row) {
                row.miss <- is.na(row)
                row.exp <- which(row.miss)
                d <- similarityCal(row[-row.exp], xcomplete[, 
                  -row.exp, drop = F], sim.method)
                id.idx <- order(d, decreasing = T)[2:(K + 1)]
                id.sel <- d[id.idx]
                const <- 1/sum(id.sel)
                w <- const * id.sel
                row[row.exp] <- w \%*\% xcomplete[id.idx, row.exp, 
                  drop = F]
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
