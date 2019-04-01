\name{fit_grp}
\alias{fit_grp}
\alias{grplasso-package}
\title{Function to fit a group lasso model to a standardized feature matrix}
\description{Standardizes feature matrix including categorical features and fits a group lasso model}
\usage{
fit_grp(eqn, dat, lambda, model = LinReg(), nonpen = c(), standardize = TRUE, ...)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{eqn}{
\code{formula} of the penalized variables. The response
    has to be on the left hand side of \code{~}. If interaction terms are included without main effects, the main effects will automatically be added by the package. 
}
  \item{dat}{\code{data.frame}, categorical features need to be of type \code{factor}
}
  \item{lambda}{
    Penalty parameter (scalar)
}
  \item{model}{
    an object of class \code{grpl.model} as defined in the package
    \code{grplasso}. 
}
  \item{nonpen}{
     \code{formula} of the nonpenalized features
}
  \item{standardize}{
   logical. If true, the design matrix of the continuous features will be centered and standardized to unit norm
}
  \item{...}{
    additional arguments to be passed to the \code{grplasso} function in
    the package of the same name.
}
}
\details{
Design matrices of the categorical features and interactions between
categorical features are centered and standardized by column-wise
scaling. After fitting a group lasso model to the standardized desgin
matrix, coefficients are re-scaled and centered to the original scale of the data. Interactions between categorical and continuous features are standardized by a singular value decomposition. 


}
\value{
  A dataframe containing the coefficients of the fitted group lasso model that have been re-scaled to the original scale of the data is returned. Coefficients of interaction terms for which no observations are included in dat are returned as NA.
}
\references{
Detmer, Felicitas J., and Martin Slawski. "A Note on Coding and Standardization of Categorical Variables in (Sparse) Group Lasso Regression." arXiv preprint arXiv:1805.06915 (2018).
}
\author{
  Felicitas Detmer, \email{fdetmer@gmu.edu}
}

\examples{
data(dattest)

#---set datatype of categorical features to factor=----
dattest$X1cut=as.factor(dattest$X1cut)
dattest$X2cut=as.factor(dattest$X2cut)
dattest$X3cut=as.factor(dattest$X3cut)

table(dattest[,c("X1cut", "X2cut", "X3cut")])

#--fit group lasso models
coefs1=fit_grp(y~X1cut * X2cut +X1cut * X3cut +X2cut * X3cut, dattest, lambda=0.5, model=LinReg())
coefs2=fit_grp(y~X1cut * X2cut +X1cut * X3cut +X2cut * X3cut, dattest, lambda=0.5, model=LinReg(),
               nonpen=~X1cut)
}
\keyword{models}
\keyword{regression}
