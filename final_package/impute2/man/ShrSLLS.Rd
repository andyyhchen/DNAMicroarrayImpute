\name{ShrSLLS}
\alias{ShrSLLS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
%%  ~~function to do ... ~~
ShrSLLS impute
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
ShrSLLS(xmiss, K = 15, sim.method = "EuDist")
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
    x.complete <- xmiss[-miss.row, ]
    x.incomplete <- xmiss[miss.row, ]
    miss.rate <- rowSums(is.na(x.incomplete))
    miss.rate.order <- order(miss.rate)
    therhold <- sum(miss.rate)/length(miss.rate)
    x.incomplete <- x.incomplete[miss.rate.order, ]
    x.imputed <- t(apply(x.incomplete, 1, function(row) {
        row.miss <- which(is.na(row))
        sim <- similarityCal(row[-row.miss], x.complete[, -row.miss], 
            sim.method)
        sim.id <- order(sim, decreasing = T)[1:K]
        x.tmp <- ginv(t(x.complete[sim.id, -row.miss, drop = FALSE])) \%*\% 
            row[-row.miss, drop = FALSE]
        x.Shrink <- as.vector((1 - (K - 2) * var(x.tmp)/((length(row) - 
            sum(is.na(row))) * sum(x.tmp^2)))) * x.tmp
        row[row.miss] <- t(x.complete[sim.id, row.miss, drop = FALSE]) \%*\% 
            x.Shrink
        if (length(row.miss) <= therhold) {
            x.complete <<- rbind(x.complete, row)
        }
        return(row)
    }))
    xmiss[miss.row[miss.rate.order], ] <- x.imputed
    return(xmiss)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
