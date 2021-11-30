\name{getCandidateEdges.WS}
\title{getCandidateEdges.WS}
\alias{getCandidateEdges}
\alias{getCandidateEdges.WS}
\description{Queries the PrixFixe web service to find all edges between
  the specified NCBI Gene IDs.
  The reference network is that used and described as part of the
  accompanying reference.
}
\usage{getCandidateEdges.WS(ncbi_gene_ids)}
\arguments{
  \item{ncbi_gene_ids}{
    A character vector of NCBI Gene IDs.
  }
  \item{edge_source}{
  a character specifying the source that was used for obtaining the
  genes that are functionally associated.
  Default value: 'MTasan_example' (The associated genes for the example_dbsnp_ids)

  For 'MTasan_paper' value it will return the gene associations Rdata described in M. Tasan paper.
  This is a data frame with almost 10 million rows.

  For 'STRING' value the function will call and parse the results from https://string-db.org/
  In this case the species is homo sapiens (9606) and STRING required_score=400.
  See also (https://string-db.org/cgi/help.pl?subpage=api%23getting-the-string-network-interactions)
  }
}
\details{
  Might return a large data frame or perform a web service call on a public server (https://string-db.org/), so
  please use judiciously/kindly!
}
\value{
  A two-column data frame holding edges connecting any of the specified
  genes within the network used by the accompanying manuscript.
  Please refer to the manuscript for greater details on the construction
  of this network, as well as the interpretation of the edges
  (i.e. ``functional associations'').
}
\references{
M Tasan, G Musso, T Hao, M Vidal, CA MacRae & FP Roth. Selecting causal genes from genome-wide association studies via functionally coherent subnetworks. Nature Methods 12(2):154-159 (2015)
}
\examples{
x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
y <- getGeneRegionsTable(x$ld_regions)
z <- getCandidateEdges.WS(y$id)
z
}
