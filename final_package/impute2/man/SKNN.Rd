\name{SKNN}
\alias{SKNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	SKNN imputation
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
SKNN(xmiss, K = 15, sim.method = "EuDist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmiss}{
%%     ~~Describe \code{xmiss} here~~
}
  \item{K}{
%%     ~~Describe \code{K} here~~
}
  \item{sim.method}{
%%     ~~Describe \code{sim.method} here~~
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
function (xmiss, K = 15, sim.method = "EuDist") 
{
    miss.gene <- is.na(xmiss)
    miss.row <- which(rowSums(miss.gene) != 0)
    xincomplete <- xmiss[miss.row, ]
    miss.inc <- is.na(xincomplete)
    miss.origin <- order(rowSums(miss.inc))
    xincomplete <- xincomplete[order(rowSums(miss.inc)), ]
    xcomplete <- xmiss[-miss.row, ]
    xtmp <- matrix(nc = ncol(xincomplete))
    xmiss[miss.row[miss.origin], ] <- t(apply(xincomplete, 1, 
        function(row) {
            row.miss <- is.na(row)
            row.exp <- which(row.miss)
            d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, 
                drop = F], sim.method)
            id.idx <- order(d, decreasing = T)[1:K]
            id.sel <- d[id.idx]
            const <- sum(id.sel)
            w <- 1/const * id.sel
            row[row.exp] <- w \%*\% xcomplete[id.idx, row.exp, 
                drop = F]
            xcomplete <<- rbind(xcomplete, row)
            return(row)
        }))
    return(xmiss)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
