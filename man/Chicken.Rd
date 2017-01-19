\name{Chicken.rda}
\alias{Chicken.rda}
\docType{data}
\title{
Sample data to run BEFDR
}
\description{
Example data needed to run function BEFDR function

}
\usage{data("Chicken")}
\format{
  A .RData object contains several vecotrs and data.frames
  \describe{
    \item{\code{sigMrk}}{ a vector of SNPs that captures the polygenic backgroud, so I would like to keep them as fixed effect, which means they will be kept in the minimal model}
    \item{\code{mrksT}}{a  set of markers I want to do BE on}
    \item{\code{mrk_geno}}{ is a genotype matrix}
    \item{\code{phe}}{ a vector contains phenotype}
    \item{\code{phe}}{ fx1 is a vector contains sex}
    \item{\code{mrk_info}}{ a matrix contains the name and location of all the marker}


  }
}
\details{
%%  ~~ If necessary, more details than the __description__ above ~~
}
\source{
A complex multi-locus, multi-allelic genetic architecture underlying the long-term selection-response in the Virginia body weight line of chickens (In prepartion)
}
\references{
Sheng Z, Pettersson ME, Honaker CF, Siegel PB, Carlborg Ö. 2015. Standing genetic variation as a major contributor to adaptation in the Virginia chicken lines selection experiment. Genome Biol 16: 219.

A complex multi-locus, multi-allelic genetic architecture underlying the long-term selection-response in the Virginia body weight line of chickens (In prepartion)

}
\examples{
data(Chicken.RData)
}
\keyword{datasets}
