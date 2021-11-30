initGA <- function()
  {
    #require(plyr)
    #require(igraph)


    PUBLIC <- list()
    PRIVATE <- list() ## just for bookkeeping. could be without PRIVATE, but makes everything more consistent to read.

    ## example raw_fitness_function(s):
    PUBLIC$edgeCountFitnessFunction <- function(vertex_sequence, graph) ecount(induced.subgraph(graph, vertex_sequence))

    ## example selection_transformation(s):
    PUBLIC$cubicTransformation <- function(x, coefficient = 1) coefficient * x^3
    PUBLIC$identityTransformation <- function(x) I(x)

    PRIVATE$mutate <- function(solution_matrix, vertex_sequence_sets, mutation_rate)
      {
        stopifnot(ncol(solution_matrix) == length(vertex_sequence_sets))
        ## loop through each vertex set (GA-gene).
        for(vertex_sequence_set_i in 1:length(vertex_sequence_sets)) {
          ## short-circuit this loop for any vertex index set with only a single vertex.
          if(length(vertex_sequence_sets[[vertex_sequence_set_i]]) == 1) next
          ## which solutions (GA-indivduals) are to be mutated?
          mutation_ix <- which(sample(c(TRUE, FALSE), size = nrow(solution_matrix), replace = TRUE, prob = c(mutation_rate, 1 - mutation_rate)))
          ## sample the new vertices (GA-alleles).
          new_vertex_indices <- vertex_sequence_sets[[vertex_sequence_set_i]][sample.int(length(vertex_sequence_sets[[vertex_sequence_set_i]]), length(mutation_ix), replace = TRUE)]
          ## update the solutions for that vertex set.
          solution_matrix[mutation_ix, vertex_sequence_set_i] <- new_vertex_indices
        }
        return(solution_matrix)
      } ## ends PRIVATE$mutate(...)


    ##------------------------------------------------------------------------------------------
    ## "fitnesses" should be directly proportional to the sampling probabilities for successful mating.
    ## this function is very slow and should almost certainly be eschewed in favor of mate.2 (see below).
    ##------------------------------------------------------------------------------------------
    PRIVATE$mate.1 <- function(solution_matrix, fitnesses)
      {
        stopifnot(nrow(solution_matrix) == length(fitnesses))
        RV <- matrix(-1, nrow = nrow(solution_matrix), ncol = ncol(solution_matrix), dimnames = dimnames(solution_matrix))
        for(mating_i in 1:nrow(solution_matrix)) {
          parent_ix <- sample.int(nrow(solution_matrix), 2, replace = FALSE, prob = fitnesses)
          RV[mating_i,] <- solution_matrix[parent_ix[1],]
          allele_swap_lix <- sample(c(TRUE, FALSE), ncol(solution_matrix), replace = TRUE)
          RV[mating_i, allele_swap_lix] <- solution_matrix[parent_ix[2], allele_swap_lix]
        }
        return(RV)
      }


    ##------------------------------------------------------------------------------------------
    ## much faster alternative is to sample 2 * N_INDIVIDUALS once, with replacement.
    ## that is, create the indices for the parents in one sampling call.
    ## the downside is that it is now possible to get a sequence of indices such that (..., x, x, ...) might occur, where both parental genotypes are identical, leading to a clonal child.
    ## i call these "clonal matings".
    ## note that a clonal mating necessarily involves identical genotypes, but identical genotype matings may also occur even in non-clonal instances.
    ## as the population diversity decreases, the clonal mating chance increases, but it may be a manageable risk, and the speed increases certainly seem to be of large-enough benefit to consider just accepting this downside.
    ##
    ## speed comparison for a solution matrix with 10e3 rows:
    ##
    ## > system.time(foo <- mate.1(solution_matrix, runif(10e3)))
    ##    user  system elapsed
    ##  14.680   0.000  14.676
    ## > system.time(foo <- mate.2(solution_matrix, runif(10e3)))
    ## [1] "2 clonal matings."
    ##    user  system elapsed
    ##   0.050   0.000   0.046
    ##Warning message:
    ##In sample.int(nrow(solution_matrix), 2 * nrow(solution_matrix),  :
    ##  Walker's alias method used: results are different from R < 2.2.0
    ##------------------------------------------------------------------------------------------
    PRIVATE$mate.2 <- function(solution_matrix, fitnesses, VERBOSE = FALSE)
      {
        ## this line was added to handle the case when all solutions are equally bad, which led to situations where the selection fitnesses were 0 for all solutions... causing the weighted sampling call (i.e. sample(..., prob=x)) to fail.
        if(all(fitnesses == 0)) fitnesses <- rep(.Machine$double.eps, length(fitnesses))

        stopifnot(nrow(solution_matrix) == length(fitnesses))
        parent_ix_matrix <- matrix(sample.int(nrow(solution_matrix), 2 * nrow(solution_matrix), replace = TRUE, prob = fitnesses), ncol = 2)
        if(VERBOSE) print(sprintf("%d clonal matings.", sum(parent_ix_matrix[,1] == parent_ix_matrix[,2])))
        allele_swap_matrix <- sample(c(TRUE, FALSE), nrow(solution_matrix) * ncol(solution_matrix), replace = TRUE)
        parent_1_solution_matrix <- solution_matrix[parent_ix_matrix[,1], , drop = FALSE]
        parent_2_solution_matrix <- solution_matrix[parent_ix_matrix[,2], , drop = FALSE]
        parent_1_solution_matrix[allele_swap_matrix] <- parent_2_solution_matrix[allele_swap_matrix]
        return(parent_1_solution_matrix)
      }


    PRIVATE$mate <- PRIVATE$mate.2


    ##------------------------------------------------------------------------------------------
    ## NOT PARALLELIZED (so can run multiple iterations in parallel for randomization cases).
    ##
    ## raw_fitness_function should take an igraph instance and return a scalar value.

    ## selection_transformation should take a numeric vector and return a numeric vector (i.e. a 1-D to 1-D transformation).
    ## stopping_statistic should also be a 1-D to a single scalar collapsing function (e.g. mean or median).
    ##------------------------------------------------------------------------------------------
    PUBLIC$run <- function(PF, n_individuals = 5e3, mutation_rate = 0.05,
                           raw_fitness_function = GPF$GA$edgeCountFitnessFunction,
                           selection_transformation = GPF$GA$cubicTransformation,
                           stopping_statistic = mean, stopping_statistic_improvement_threshold = 0.01,
                           max_n_generations = 50,
                           VERBOSE = FALSE)
      {
        n_individuals <- as.integer(n_individuals)
        stopifnot(n_individuals > 1)

        mutation_rate <- as.numeric(mutation_rate)
        stopifnot(0 <= mutation_rate && mutation_rate <= 1)


        computeFitnesses <- function(solution_matrix, G) apply(solution_matrix, 1, raw_fitness_function, G)


        ##------------------------------------------------------------------------------------------
        ## initialize the GA-population matrix.
        ## each row represents a GA-individual (GA-genotype), corresponding to a single solution.
        ## each col represents a GA-locus, corresponding to one of the simplifed gene sets.
        ##------------------------------------------------------------------------------------------
        solution_matrix <- sapply(PF$SIMPLIFIED_VERTEX_SEQUENCE_SETS, function(vertex_sequence) {
          vertex_sequence[sample.int(length(vertex_sequence), n_individuals, replace = TRUE)]
        })


        ##------------------------------------------------------------------------------------------
        ## the actual genetic algorithm.
        ## really only consists of two steps: mutate and mate (the mating encapsulates the concept of selection, via preferential matings for fitter individuals).
        ##------------------------------------------------------------------------------------------
        fitnesses_per_generation <- list()
        n_distinct_solutions_per_generation <- list()

        fitnesses_per_generation[[1]] <- computeFitnesses(solution_matrix, PF$G)
        n_distinct_solutions_per_generation[[1]] <- nrow(unique(solution_matrix))

        if(VERBOSE) print("generation 0:")
        if(VERBOSE) print(summary(fitnesses_per_generation[[0 + 1]]))

        for(foo.generation_i in 1:max_n_generations) {
          solution_matrix <- PRIVATE$mutate(solution_matrix, PF$SIMPLIFIED_VERTEX_SEQUENCE_SETS, mutation_rate)
          solution_matrix <- PRIVATE$mate(solution_matrix, selection_transformation(computeFitnesses(solution_matrix, PF$G)))

          fitnesses_per_generation[[foo.generation_i + 1]] <- computeFitnesses(solution_matrix, PF$G)
          n_distinct_solutions_per_generation[[foo.generation_i + 1]] <- nrow(unique(solution_matrix))

          if(VERBOSE) print(sprintf("generation %d:", foo.generation_i))
          if(VERBOSE) print(summary(fitnesses_per_generation[[foo.generation_i + 1]]))

          ## if M(i) is the median fitness for generation i, then the relative improvement of generation i vs generation i-1 is:
          ## (M(i) - M(i-1)) / M(i-1)
          stopping_improvement <- NA2((stopping_statistic(fitnesses_per_generation[[foo.generation_i + 1]]) - stopping_statistic(fitnesses_per_generation[[foo.generation_i]])) / stopping_statistic(fitnesses_per_generation[[foo.generation_i]]), 0)
          if(VERBOSE) print(sprintf("relative stopping-statistic improvement over previous generation = %G", stopping_improvement))
          if(stopping_improvement < stopping_statistic_improvement_threshold) break
        }


        RV <- list()
        RV$solution_matrix <- solution_matrix
        RV$fitnesses_per_generation <- fitnesses_per_generation
        RV$n_distinct_solutions_per_generation <- n_distinct_solutions_per_generation
        RV$mutation_rate <- mutation_rate
        RV$raw_fitness_function <- raw_fitness_function
        RV$PF <- PF
        class(RV) <- "GARunResult"
        return(RV)
      }


    PUBLIC$scoreVertices.retired_2013_10_09 <- function(GARunResult, return_as_table = FALSE, .parallel = FALSE, VERBOSE = FALSE)
      {
        stopifnot(is(GARunResult, "GARunResult"))

        SOLUTION_MATRIX <- GARunResult$solution_matrix
        SIMPLIFIED_VERTEX_SEQUENCE_SETS <- GARunResult$PF$SIMPLIFIED_VERTEX_SEQUENCE_SETS
        G <- GARunResult$PF$G
        FITNESS_FUNCTION <- GARunResult$fitness_function

        ## the igraph vertex sequence (scalar) for the "foo" vertex.
        foo_vertex_sequence1 <- which(V(G)$name == "foo")

        ##------------------------------------------------------------------------------------------
        ## loop through each vertex set (GA-gene).
        ## represents a single column in the solution_matrix.
        ## we'll modify solution_matrix here, so either be sure it's in a function environment (and thus local to that new environment), or is explicitly made local if a for() loop is used.
        ## basically, be sure not to screw up the actual solution_matrix object during these manipulations.
        ##------------------------------------------------------------------------------------------
        vertex_scores <- llply(1:length(SIMPLIFIED_VERTEX_SEQUENCE_SETS), function(vertex_set_i) {
          if(VERBOSE) print(sprintf("scoring vertex set %s.", names(SIMPLIFIED_VERTEX_SEQUENCE_SETS)[vertex_set_i]))

          ## score the "foo" solution for this GA-gene (i.e. force this GA-gene to have the "foo" GA-allele).
          ## (a vector of length == nrow(solution_matrix)).
          SOLUTION_MATRIX[,vertex_set_i] <- rep(foo_vertex_sequence1, nrow(SOLUTION_MATRIX))
          foo_fitnesses <- apply(SOLUTION_MATRIX, 1, FITNESS_FUNCTION, G)

          ## now score each other vertex (candidate) in turn (same forcing idea as above, just with each actual GA-allele).
          ## nrow(candidate_edge_counts) == nrow(solution_matrix).
          ## ncol(candidate_edge_counts) == length(vertex_sequence_sets[[vertex_set_i]]); i.e. a column for each GA-allele in this vertex set.
          candidate_fitnesses <- sapply(SIMPLIFIED_VERTEX_SEQUENCE_SETS[[vertex_set_i]], function(vertex_sequence1) {
            SOLUTION_MATRIX[,vertex_set_i] <- rep(vertex_sequence1, nrow(SOLUTION_MATRIX))
            apply(SOLUTION_MATRIX, 1, FITNESS_FUNCTION, G)
          })

          ## for each solution, what's the best edge count seen (keeping all other GA-genes) fixed to their GA-found solutions?
          max_fitnesses <- apply(candidate_fitnesses, 1, max)

          ## what's the difference between the max and the min?
          ## note that the min is equivalent to the "foo" vector.
          max_minus_foo_fitnesses <- max_fitnesses - foo_fitnesses

          ## now linearly rescale each candidate edge count from [foo,max] to [0,1].
          ## note the zero-variance case will introduce two possibilities, described further in the notes.
          ## the NA2(x,0) case handles the zero-variance scenario where foo == x == max.
          linearly_rescaled_candidate_fitnesses <- apply(candidate_fitnesses, 2, function(x) {
            NA2((x - foo_fitnesses) / (max_minus_foo_fitnesses), 0)
          })

          ## finally incorporate the GA-gene and solution scores, as described in the notes, and average.
          RV <- apply(linearly_rescaled_candidate_fitnesses, 2, function(x) mean(x * max_minus_foo_fitnesses))
          names(RV) <- V(G)$name[SIMPLIFIED_VERTEX_SEQUENCE_SETS[[vertex_set_i]]]
          return(RV)
        }, .parallel = .parallel)
        names(vertex_scores) <- names(SIMPLIFIED_VERTEX_SEQUENCE_SETS)

        if(!return_as_table) {
          return(vertex_scores)
        } else {
          rv <- do.call(rbind, lapply(names(vertex_scores), function(i) {
            data.frame(set_name = rep(i, length(vertex_scores[[i]])),
                       vertex_name = names(vertex_scores[[i]]),
                       score = vertex_scores[[i]],
                       stringsAsFactors = FALSE)
          }))
          rownames(rv) <- NULL
          return(rv)
        }
      }


    PUBLIC$scoreVertices <- function(GARunResult, return_as_table = TRUE, .parallel = FALSE, VERBOSE = FALSE)
      {
        stopifnot(is(GARunResult, "GARunResult"))

        SOLUTION_MATRIX <- GARunResult$solution_matrix
        SIMPLIFIED_VERTEX_SEQUENCE_SETS <- GARunResult$PF$SIMPLIFIED_VERTEX_SEQUENCE_SETS
        G <- GARunResult$PF$G
        RAW_FITNESS_FUNCTION <- GARunResult$raw_fitness_function

        ## the igraph vertex sequence (scalar) for the "foo" vertex.
        foo_vertex_sequence1 <- which(V(G)$name == "foo")

        ##------------------------------------------------------------------------------------------
        ## loop through each vertex set (GA-gene).
        ## represents a single column in the solution_matrix.
        ## we'll modify solution_matrix here, so either be sure it's in a function environment (and thus local to that new environment), or is explicitly made local if a for() loop is used.
        ## basically, be sure not to screw up the actual solution_matrix object during these manipulations.
        ##------------------------------------------------------------------------------------------
        vertex_scores <- llply(1:length(SIMPLIFIED_VERTEX_SEQUENCE_SETS), function(vertex_set_i) {

          if(VERBOSE) print(sprintf("scoring vertex set %s.", names(SIMPLIFIED_VERTEX_SEQUENCE_SETS)[vertex_set_i]))

          ## score the "foo" solution for this GA-gene (i.e. force this GA-gene to have the "foo" GA-allele).
          ## (a vector of length == nrow(solution_matrix)).
          SOLUTION_MATRIX[,vertex_set_i] <- rep(foo_vertex_sequence1, nrow(SOLUTION_MATRIX))
          foo_fitnesses <- apply(SOLUTION_MATRIX, 1, RAW_FITNESS_FUNCTION, G)

          ## now score each other vertex (candidate) in turn (same forcing idea as above, just with each actual GA-allele).
          ## nrow(candidate_fitnesses) == nrow(solution_matrix).
          ## ncol(candidate_fitnesses) == length(vertex_sequence_sets[[vertex_set_i]]); i.e. a column for each GA-allele in this vertex set.
          candidate_fitnesses <- sapply(SIMPLIFIED_VERTEX_SEQUENCE_SETS[[vertex_set_i]], function(vertex_sequence1) {
            SOLUTION_MATRIX[,vertex_set_i] <- rep(vertex_sequence1, nrow(SOLUTION_MATRIX))
            apply(SOLUTION_MATRIX, 1, RAW_FITNESS_FUNCTION, G)
          })

          ## difference from foo vertex.
          difference_score <- candidate_fitnesses - foo_fitnesses

          RV <- apply(difference_score, 2, mean)
          names(RV) <- V(G)$name[SIMPLIFIED_VERTEX_SEQUENCE_SETS[[vertex_set_i]]]
          return(RV)
        }, .parallel = .parallel)
        names(vertex_scores) <- names(SIMPLIFIED_VERTEX_SEQUENCE_SETS)

        ##------------------------------------------------------------------------------------------
        ## regardless of the chosen solution's fitness, the maximum difference beteen the foo vertex and a candidate vertex is length(SIMPLIFIED_VERTEX_SEQUENCE_SETS) - 1.
        ## this maximum case only occurs if the candidate vertex is connected to every other set.
        ## we can scale by this maximum to get vertex scores that are independent of the graph size.
        ## this will still preserve solution importance, since the lower-scoring solutions have fewer edges period, the difference between the candidates and foo should also be smaller.
        ##------------------------------------------------------------------------------------------
        scaled_vertex_scores <- llply(vertex_scores, function(x) x / (length(SIMPLIFIED_VERTEX_SEQUENCE_SETS) - 1))


        if(!return_as_table) {
          return(vertex_scores)
        } else {
          rv <- do.call(rbind, lapply(names(vertex_scores), function(i) {
            data.frame(set_name = rep(i, length(vertex_scores[[i]])),
                       vertex_name = names(vertex_scores[[i]]),
                       score = vertex_scores[[i]],
                       scaled_score = scaled_vertex_scores[[i]],
                       stringsAsFactors = FALSE)
          }))
          rownames(rv) <- NULL
          return(rv)
        }

      }

    return(PUBLIC)
  }

