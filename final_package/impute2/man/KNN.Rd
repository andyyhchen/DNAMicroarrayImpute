\name{KNN}
\alias{KNN}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	KNN DNA impute 
}

\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
KNN(x, K = 15, sim.method = "EuDist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{x}{
	  matirx with missing value
}
  \item{K}{
	  number of selected neighbor 
}
  \item{sim.method}{
	  method of neighbors selection 
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
function (x, K = 15, sim.method = "EuDist") 
{
    miss.gene <- is.na(x)
    miss.row <- which(rowSums(miss.gene) != 0)
    xcomplete <- x[-miss.row, ]
    xincomplete <- x[miss.row, ]
    x[miss.row, ] <- t(apply(xincomplete, 1, function(row) {
        row.miss <- is.na(row)
        row.exp <- which(row.miss)
        d <- similarityCal(row[-row.exp], xcomplete[, -row.exp, 
            drop = F], sim.method)
        id.idx <- order(d, decreasing = T)[1:K]
        id.sel <- d[id.idx]
        const <- sum(id.sel)
        w <- 1/const * id.sel
        w <- matrix(w, nc = length(w), nr = 1)
        row[row.exp] <- w \%*\% xcomplete[id.idx, row.exp, drop = F]
        return(row)
    }))
    return(x)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
