\name{LS}
\alias{LS}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
	LS impute
}
\description{
	LS DNA imputation 
}
\usage{
LS(xmiss, K = 15, sim.method = "cor")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{xmiss}{
	  matrix with missing value
}
  \item{K}{
	  numbers of neighbors 
}
  \item{sim.method}{
	  method of choosing neighbors
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
function (xmiss, K = 15, sim.method = "cor") 
{
    miss.gene <- is.na(xmiss)
    miss.row <- which(rowSums(miss.gene) != 0)
    xcomplete <- xmiss[-miss.row, ]
    xincomplete <- xmiss[miss.row, ]
    impute <- function(row) {
        row.exp <- which(is.na(row))
        gene <- row[-row.exp]
        cand_x <- xcomplete[, -row.exp, drop = F]
        sim <- similarityCal(gene, cand_x, method = sim.method)
        row.idx <- order(sim, decreasing = T)[1:K]
        row.r <- sim[row.idx]
        row.cand <- cand_x[row.idx, , drop = F]
        lg <- apply(row.cand, 1, function(x) {
            lm(gene ~ x)$coefficients
        })
        row.impcand <- xcomplete[row.idx, row.exp, drop = F]
        y <- matrix(0, nc = ncol(row.impcand), nr = nrow(row.impcand))
        w <- (row.r^2/(1 - row.r^2 + 1e-06))^2
        sw <- sum(w)
        w <- w/sw
        for (i in 1:nrow(row.impcand)) {
            y[i, ] <- lg[2, i] * row.impcand[i, ] + lg[1, i]
        }
        row[row.exp] <- apply(y, 2, function(x) {
            sum(w * x)
        })
        return(row)
    }
    xmiss[miss.row, ] <- t(apply(xincomplete, 1, impute))
    return(xmiss)
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
