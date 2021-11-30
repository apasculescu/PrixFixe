\name{GPF$GA$run}
\title{GPF$GA$run}
\alias{GPF$GA$run}
\description{Runs the genetic algorithm optimizer to find dense prix
  fixe subnetworks.}
\usage{GPF$GA$run(PF,
  n_individuals = 5e3,
  mutation_rate = 0.05,
  raw_fitness_function = GPF$GA$edgeCountFitnessFunction,
  selection_transformation = GPF$GA$cubicTransformation,
  stopping_statistic = mean, stopping_statistic_improvement_threshold = 0.01,
  max_n_generations = 50,
  VERBOSE = FALSE)
}
\arguments{
  \item{PF}{
    A prix fixe object, as constructed by the GPF$PF$new function.
  }
  \item{n_individuals}{
    The number of individuals in the genetic algorithm's population.
    (default = 5000).
  }
  \item{mutation_rate}{
    The mutation probability/generation (between 0 and 1) of each individual in the genetic algorithm
    population (default = 0.05).
  }
  \item{raw_fitness_function}{
    A fitness function by which to judge the current candidate solution
    (i.e. "individual").
    Currently, the only supported fitness function is edge count
    (i.e. unweighted density), which is set as the default.
    It is not recommended to try to alter this argument just yet, as
    the function's specifications are subject to change.
  }
  \item{selection_transformation}{
    A vectorized transformation function that takes a numeric vector and
    returns a numeric vector.
    The returned vector is used to help "upweight" fit individuals
    during the genetic algorithm's selection process.
    Details are available in the accompanying manuscript, and the
    default value is a cubic transformation (also described in the
    manuscript)
  }
  \item{stopping_statistic}{
    A reducing function that reports a single statistic for the
    population at each generation by which that generation's aggregate
    fitness is assessed.
    The default is the mean fitness.
  }
  \item{stopping statistic_improvement_threshold}{
    The relative incremental value of the stopping_statistic at which
    the algorithm will terminate.
    That is, when the relative improvement over the previous generation
    is below this threshold, the collective fitness is said to have
    "plateau-ed".
    The default value is 0.01.
  }
  \item{max_n_generations}{
    In case the fitness doesn't plateau in reasonable time, this is an
    upper-bound on the number of generations used in the genetic
    algorithm.
    The default is 50, but this default is rarely reached except by the
    largest of realistic inputs.
    And if you're inputing that many loci, we should probably talk to
    see what other optimizations might be possible.
  }
  \item{VERBOSE}{
    A logical value (with FALSE as the default) that allows for
    increased reporting as the function is running (e.g. tracking
    fitness statistics during generational progression).
  }
}
\details{
  This function runs the genetic algorithm to find dense prix fixe
  subnetworks.
  This is a CPU-intensive task, and may take a while to complete.
  Use VERBOSE = TRUE to follow along while it runs (and to convince
  yourself it hasn't frozen).
}
\value{
  An object of class "GARunResult", which is a list with many
  componenets, the following of which are of most interest:
  \itemize{
    \item{PF}{
      A copy of the PrixFixe object passed as the first argument.
      This copy is needed for internal referencing in case the original
      PrixFixe object is modified, since the algorithm's run results use
      vertex sequences for the underying graph object, which is stored
      in the PrixFixe object.
    }
    \item{solution_matrix}{
      This is a matrix with a row for each individual in the final
      population, and a column for each locus.
      The value in a cell corresponds to the vertex sequence identifier
      in the PF$G (igraph) object for that individual.
      This, each row represents a single prix fixe subgraph.
    }
    \item{fitnesses_per_generation}{
      A list of vectors, where each vector is of length equal to the
      population size.
      Each vector itself represents a generation.
      The i-th value of the j-th vector represents the fitness of the
      i-th individual during the j-th generation.
    }
    \item{n_distinct_solutions_per_generation}{
      A list of values, one for each generation, that measues the
      "complexity" of the population.
      For some small inputs, the number of possible solutions dictates
      that many individuals during the algorithm's progression will be
      identical with respect to each other.
      The mutation rate and population size can be altered to force
      higher (or lower) complexity, if needed.
    }
  }
  This object isn't intended for by-hand processing, but may be
  interesting to examine.
  To obtain gene scores from this object, use the GPF$GA$scoreVertices function.
}
\examples{
  x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
  y <- getGeneRegionsTable(x$ld_regions)
  z <- getCandidateEdges.WS(y$id)
  pf <- GPF$PF$new(tapply(y$id, y$region_name, I), z, VERBOSE = TRUE)
  run_result <- GPF$GA$run(pf, VERBOSE = TRUE)
}

