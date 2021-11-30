\name{getGeneRegionsTable}
\title{getGeneRegionsTable}
\alias{getGeneRegionsTable}
\description{Extracts gene data from an ld_region object into a single
  table.}
\usage{getGeneRegionsTable(ld_regions)}
\arguments{
  \item{ld_regions}{
    An ld_regions object, which is typically extracted from a call to
    getLDRegionsFromSNPs.WS.
    See example below.
  }
}
\details{
  This function is just a convenience method to simplfy the results from
  an ld_regions object.
  Often people would rather see all candidate genes across all loci in a
  single table, and this function creates that table.
}
\value{
  A data frame with all genes across all loci represented in the
  specified ld_regions object.
}
\examples{
x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
y <- getGeneRegionsTable(x$ld_regions)
y
}
