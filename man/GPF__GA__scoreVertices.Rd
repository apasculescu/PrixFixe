\name{GPF$GA$scoreVertices}
\title{GPF$GA$scoreVertices}
\alias{GPF$GA$scoreVertices}
\description{Scores genes after having run the genetic algorithm prix
  fixe optimizer.}
\usage{GPF$GA$scoreVertices(GARunResult,
  return_as_table = TRUE,
  .parallel = FALSE,
  VERBOSE = FALSE)
}
\arguments{
  \item{GARunResult}{
    A GARunResult object, as returned by the GPF$GA$run function.
  }
  \item{return_as_table}{
    Return the results in a table (the default) or a list?
  }
  \item{.parallel}{
    Parallelize the scoring?
    This feature isn't fully functional yet, and the user can experiment
    with turning it on, but it's not guaranteed that parallelization
    support will work on the user's system.
  }
  \item{VERBOSE}{
    A logical value (with FALSE as the default) that allows for
    increased reporting as the function is running.
  }
}
\details{
  Scores the genes after having run the genetic algorithm prix fixe
  optimizer.
  See the Accompanying manuscript for details on this scoring procedure.
}
\value{
  If return_as_table is TRUE, then the result is a data frame with the
  gene's scores.
  If return_as_table is FALSE, then the result is a list of scores, one
  vector for each locus.
  IMPORTANT(!): not all of the genes passed to the GPF$PF$new
  constructor will be present here.
  Only those genes with potentially non-zero scores are represented
  here.
  Those genes with scores guaranteed to be zero have already been pruned
  during the GPF$PF$new construction process.
}
\examples{
  x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
  y <- getGeneRegionsTable(x$ld_regions)
  z <- getCandidateEdges.WS(y$id)
  pf <- GPF$PF$new(tapply(y$id, y$region_name, I), z, VERBOSE = TRUE)
  run_result <- GPF$GA$run(pf, VERBOSE = TRUE)
  gene_scores <- GPF$GA$scoreVertices(run_result, VERBOSE = TRUE)
  gene_scores
}

