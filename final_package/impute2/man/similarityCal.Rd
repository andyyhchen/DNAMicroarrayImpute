\name{similarityCal}
\alias{similarityCal}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
calculate similarity
%%  ~~function to do ... ~~
}
\description{
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
similarityCal(vec, mat, method = "EuDist")
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{vec}{
%%     ~~Describe \code{vec} here~~
}
  \item{mat}{
%%     ~~Describe \code{mat} here~~
}
  \item{method}{
%%     ~~Describe \code{method} here~~
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
function (vec, mat, method = "EuDist") 
{
    methods <- c("EuDist", "cor", "Angle")
    switch(match.arg(method, methods), EuDist = 1/sqrt(rowSums((mat - 
        matrix(vec, nc = length(vec), nr = nrow(mat), byrow = TRUE))^2)), 
        cor = apply(mat, 1, function(i) abs(cor(i, vec, use = "everything", 
            method = "pearson"))), Angle = apply(mat, 1, function(i) abs(sum(i * 
            vec)/sqrt(sum(i^2) * sum(vec^2)))))
  }
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
