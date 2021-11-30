\name{GPF$PF$new}
\title{GPF$PF$new}
\alias{GPF$PF$new}
\description{Create a new prix fixe instance.}
\usage{GPF$PF$new(candidate_gene_sets,
  candidate_edges,
  FAN_vertices = NULL,
  VERBOSE = FALSE)}
\arguments{
  \item{candidate_gene_sets}{
    A list of character vectors of gene IDs.
    Each vector corresponds to a separate locus.
    Preferrably, the list should be named (with locus) names, but if not
    automatic names will be assigned to each locus.
  }
  \item{candidate_edges}{
    A two-column data frame or matrix containing the gene IDs of edges
    within the functional association network that connect the candidate
    genes.
    Technically, a full network can be passed here, but it's much more
    efficient to start with a network that only contains edges between
    those candidate genes list in the candidate_gene_sets argument.
    These edges will further be pruned to create the n-partite graph.
  }
  \item{FAN_vertices}{
    Unused in this release. Support for this will be added soon. This
    paramter helps to differentiate between genes that have a score of 0
    due to lack of connectivity, and those for which scores cannot be
    assigned (due to lack of that gene's appearance in the complete
    network).
  }
  \item{VERBOSE}{
    A logical value (with FALSE as the default) that allows for
    increased reporting as the funtion is running.
  }
}
\details{
  This function creates a new prix fixe instance.
  A prix fixe instance is not a solved problem, but is instead a problem
  waiting to be solved.
  The constraints are established (for example by creating a n-partite
  graph in which to search), and some early deterministic optimizations
  are performed to reduce the search space.
  (For example, loci with genes that don't connect to any other
  candidate genes are pruned from the search space.)
}
\value{
  An object of class "PrixFixe", which is a list with the following
  components:
  \itemize{
    \item{CANDIDATE_TABLE}{
      An IRanges::DataFrame object holding information about all the
      candidate genes.
    }
    \item{G}{
      An igraph graph object representing the pruned n-partite graph in
      which the dense prix fixe subgraph search will take place.
    }
    \item{SIMPLIFIED_VERTEX_SEQUENCE_SETS}{
      A list of vertex sequences (using igraph terminology) that
      specifies the source locus ("set") for each gene ("vertex").
    }
  }
  The prix fixe object isn't intended for by-hand processing, but may be
  interesting to examine (particularly the CANDIDATE_TABLE).
}
\examples{
  x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
  y <- getGeneRegionsTable(x$ld_regions)
  z <- getCandidateEdges.WS(y$id)
  pf <- GPF$PF$new(tapply(y$id, y$region_name, I), z, VERBOSE = TRUE)
  str(pf)
}

