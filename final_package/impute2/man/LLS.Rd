\name{LLS}
\alias{LLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
LLS impute
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
LLS(xmiss, K = 15, sim.method = "EuDist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmiss}{
	  matrix with missing value
}
  \item{K}{
	  numbers of neighbors choosen 
}
  \item{sim.method}{
	  method of selectng neighbors
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
    x.complete <- xmiss[-miss.row, ]
    x.incomplete <- xmiss[miss.row, ]
    x.imputed <- t(apply(x.incomplete, 1, function(row) {
        row.miss <- which(is.na(row))
        sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], 
            sim.method)
        sim.id <- order(sim, decreasing = T)[1:K]
        row[row.miss] <- t(x.complete[sim.id, row.miss, drop = FALSE]) \%*\% 
            ginv(t(x.complete[sim.id, -row.miss, drop = FALSE])) \%*\% 
            row[-row.miss, drop = FALSE]
        return(row)
    }))
    xmiss[miss.row, ] <- x.imputed
    return(xmiss)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
