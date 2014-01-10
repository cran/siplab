\encoding{UTF-8}
\name{siplab-package}
\alias{siplab-package}
\alias{siplab}
\docType{package}
\title{
Spatial Individual-Plant Simulation
}
\description{
A platform for experimenting with spatially explicit individual-based plant modelling
}
\details{
\tabular{ll}{
Package: \tab siplab\cr
Type: \tab Package\cr
Version: \tab 1.0\cr
Date: \tab 2014-01-09\cr
License: \tab GPL\cr
}
%~~ An overview of how to use the package, including the most important functions ~~
The main function is \code{\link{assimilation}}, which computes effective resource capture (or \dQuote{assimilation}) indices. One starts with a spatial resource distribution that is typically assumed to be uniform, Plants exert competitive pressure depending on size and distance, described by \link{influence} functions. The resource available at each point is allocated to plants according to their local influence and to a partition rule. Finally, the resource uptake is weighted by an \link{efficiency} function that depends on size and distance, and is spatially integrated to obtain the plant's assimilation index.

Several examples of influence and efficiency functions are pre-programmed, see \link{influence} and \link{efficiency}. Some sample data sets are included, see links below.

Another class of spatial individual-plant models consider only pair-wise interactions between a plant and each of its neighbors, ignoring higher-order interactions. These are implemented in function \code{\link{pairwise}}.

The \code{\link{edges}} function is useful for handling edge effects.

The package is built on top of the \pkg{spatstat} library (\url{http://www.spatstat.org}), which needs to be installed first.
}
\author{
Oscar \enc{García}{Garcia}

Maintainer: O. Garcia <garcia@unbc.ca>
%~~ The author and/or maintainer of the package ~~
}
\references{
%~~ Literature or other references for background information ~~
\url{http://forestgrowth.unbc.ca/siplab}

\enc{García}{Garcia}, O. \dQuote{A generic approach to spatial individual-based modelling and simulation of plant communities}. Manuscript in review.
}
%~~ Optionally other standard keywords, one per line, from file KEYWORDS in the R documentation directory ~~
%\keyword{ package }
\seealso{
%~~ Optional links to other man pages, e.g. ~~
%~~ \code{\link[<pkg>:<pkg>-package]{<pkg>}} ~~
Example \pkg{siplab} data sets: \code{\link{boreasNP}}, \code{\link{boreasNS}}, \code{\link{boreasSA}}, \code{\link{boreasSP}}.

Some \pkg{spatstat} standard data sets may also be of interest: \code{finpines}, \code{longleaf}, \code{spruces}, \code{waka}.
}
\examples{
%~~ simple examples of the most important functions ~~
a <- assimilation(finpines, influence=zoi.inf, infpar=c(k=0.2,
       smark=1), partpar=1)
head(marks(a))
}

\section{Warning}{
Limited testing so far, use with care. Please email issues to \email{garcia@unbc.ca}.
}
