\name{gen.copula}
\alias{gen.copula}

\title{A function to generate appropriate functions to be used to fit a bivariate copula in GAMLSS}
\description{
	The \code{gen.copula()} function allows the user to generate a \code{d}, \code{p}, \code{q}, \code{p2c1}, \code{q2c1} , \code{plot} and fitting \code{gamlss} functions for bivariate copulas.
	The function can take any \code{\link[gamlss.dist]{gamlss.family}} distribution.
}
\usage{
	gen.copula(m1 = "NO", m2 = "NO", copula = "NO",
		name = NULL, print = TRUE)
}

\arguments{
	\item{m1}{character, a \code{\link[gamlss.dist]{gamlss.family}} distribution, which is 
		used to define the first marginal distribution.
		The distribution families supported by \code{gamlss2()} can be found in \code{\link[gamlss.dist]{gamlss.family}} and in the package \code{gamlss.dist}.}
	\item{m2}{character, a \code{\link[gamlss.dist]{gamlss.family}} distribution, which is 
		used to define the second marginal distribution.}
	\item{copula}{character, a copula family. For currently supported families, see details.}
	\item{name}{character, the name of the new functions. The default is a concatenation of \code{m1}, \code{m2} and \code{copula}.}
	\item{print}{boolean, whether to print the function execution (TRUE) or not (FALSE)}
}

\value{
  Returns the \code{d} (i.e., probability density function [pdf]), \code{p} (i.e., cumulative distribution function [cdf]), \code{r} (i.e., random sample generator),
	\code{p2c1} (i.e., cdf of the second marginal given the first marginal), \code{q2c1} (i.e., inverse cdf of the second marginal given the first marginal),
	\code{plot}, and the fitting function used for \code{gamlss2()}.
}

\details{
	Copulas can have positive dependence [pd], negative dependence [nd], or both.

	The currently supported one-parameter copula families are: 
	\itemize{
		\item \code{NO}: Gaussian (normal) copula [pd/nd]
		\item \code{PL}: Plackett copula [pd/nd]
		\item \code{FR}: Frank copula [pd/nd]
		\item \code{MTCJ}: Mardia-Takahasi-Clayton-Cook-Johnson copula [pd]
		\item \code{JOE}: Joe copula [pd]
		\item \code{GU}: Gumbel copula [pd]
		\item \code{GAL}: Galambos copula [pd]
		\item \code{HR}: Husler-Reiss copula [pd]
		\item \code{FGM}: Farlie-Gumbel-Morgenstern copula [pd/nd]
	}
	The currently supported two-parameter copula families are: 
	\itemize{
		\item \code{T}: Student copula [pd/nd]
		\item \code{BB1}: BB1 copula [pd]
	}
	
	For a fitted gamlss2 model with a bivariate copula family, five types of quantile
	residuals are available and can be selected with the \code{id} parameter:
	\itemize{
		\item (\code{id = 1}) First marginal
		\item (\code{id = 2}) Second marginal
		\item (\code{id = 3}) Second marginal conditional on first marginal
		\item (\code{id = 4}) Chi-square distribution with two degrees of freedom induced by p, p = res(id=1)^2 + res(id=3)^2
		\item (\code{id = 5}, default) Joint residuals, according to Equation 5 of Kalliovirta (2008)
	}
	If \code{id = 0} is used, all the five options above will be returned in a matrix format.
}

\references{
	Kalliovirta, L. (2008). \emph{Quantile residuals for multivariate models}. Helsinki Center of
	Economic Research Discussion Papers, 247.

	Joe, H. (2023). \emph{Dependence Modeling with Copulas}. CRC Press. ISBN 978-1032477374.

	Stasinopoulos D. M., Rigby R.A., Heller G., Voudouris V., and De Bastiani F., (2017)
	\emph{Flexible Regression and Smoothing: Using GAMLSS in R},  Chapman and Hall/CRC.  
	(see also \url{https://www.gamlss.com/}).
}
\author{Arthur Fendrich \email{arthur.fendrich@ec.europa.eu}, Elise Van Eynde, Mikis Stasinopoulos and Bob Rigby}

\examples{
library(copula)
library(gamlss2)
library(gamlss.dist)
library(gamlss.copula)

## Generates two log-modified families
gen.Family('SN2', type = 'log')
gen.Family('SEP4', type = 'log')

## Generates the bivariate copula
gen.copula('logSN2', 'logSEP4', 'PL', name = 'copula')

## Defines the parameter formulas
eq.mu1 <- y ~ 1
eq.sigma1 <- eq.nu1 <- ~ 1
eq.mu2 <- eq.sigma2 <- eq.nu2 <- eq.tau2 <- ~ 1
eq.theta <- ~ 1
eq.lst <- list(eq.mu1, eq.sigma1, eq.nu1,
	eq.mu2, eq.sigma2, eq.nu2, eq.tau2,
	eq.theta)

## Simulates some data with known parameters
set.seed(1)
xy <- rCopula(2500, plackettCopula(param = 3))
y.m1 <- qlogSN2(xy[,1], mu = 1, sigma = 0.75, nu = 0.5)
y.m2 <- qlogSEP4(xy[,2], mu = 0, sigma = 1.2, nu = 2.4, tau = 1)

## Estimates the model
n <- length(y.m1)
y <- c(y.m1, y.m2)
b <- gamlss2(eq.lst, data = data.frame(y = y), family = copula,
	control = gamlss2_control(maxit = c(1000, 10)))

## Validation below
predict(b)[1,]

## Plots log-density of this copula
plot.copula(m1.mu = 1, m1.sigma = 0.75, m1.nu = 0.5,
	m2.mu = 0, m2.sigma = 1.2, m2.nu = 2.4, m2.tau = 1,
	theta = 3, log = TRUE, plot = TRUE)
}
\keyword{distribution}
\keyword{regression}
\keyword{bivariate}
\keyword{copula}
