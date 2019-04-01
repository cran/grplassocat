\name{internals}
\alias{gen_dum}
\alias{func_int}
\alias{seqWrapper}
\alias{dattest}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{Miscellaenous}
\description{\code{gen_dum}, \code{func_int}, and \code{seqWrapper} are internal functions not intended to be called directly by the
user. \code{dattest} is a synthetic data set for testing and
demonstration purposes.  
}
%\usage{
%gen_dum(dat)
%}
%- maybe also 'usage' for other objects documented here.
%\arguments{
%%     ~~Describe \code{dat} here~~
%}
%}
%\details{
%%  ~~ If necessary, more details than the description above ~~
%}
%\value{
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
%}
%\references{
%% ~put references to the literature/web site here ~
%}
\author{
Felicitas Detmer, \email{fdetmer@gmu.edu}
}
%\note{
%%  ~~further notes~~
%}

%% ~Make other sections like Warning with \section{Warning }{....} ~

%\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
%}
%\examples{
%##---- Should be DIRECTLY executable !! ----
%##-- ==>  Define data, use random,
%##--	or do  help(data=index)  for the standard data sets.
%
%## The function is currently defined as
%function (dat) 
%{
%    x = list()
%    for (i in seqWrapper(1, ncol(dat))) {
%        lvls = as.character(sort(unique(dat[, i])))
%        cn = colnames(dat)[i]
%        dums = data.frame(matrix(0, nrow = nrow(dat), ncol = length(lvls)))
%        dat[, i] = as.character(dat[, i])
%        colnames(dums) = paste(cn, lvls, sep = ".")
%        for (j in 1:nrow(dat)) {
%            dums[j, paste(cn, dat[j, i], sep = ".")] = 1
%        }
%        x[[i]] <- dums
%    }
%    if (length(x) > 0) {
%        return(x)
%    }
%    else {
%        return(data.frame()[1:nrow(dat), ])
%    }
%  }
%}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
%\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
%\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
