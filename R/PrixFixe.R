initPF <- function()
  {
#    require(igraph)
#    require(IRanges)
    stopifnot(exists("GPF"))
    stopifnot(sapply(c("ID", "UTIL"), function(x) {
      exists(x, envir = as.environment(get("GPF")), inherits = FALSE)
    }))


    PUBLIC <- list()


    ##------------------------------------------------------------------------------------------
    ## NOTE that the term "vertex sequence" refers to the igraph conotation of "sequence", which is a vector of IDs/indices in the igraph graph object.
    ##
    ## NOTE: this constructor does not enforce true N-partite-ness.
    ## it is up to the caller of this function to ensure that the candidate_gene_sets are all pairwise-disjoint.
    ## the PrixFixe problem is technically still 'valid' if a gene appears in multiple sets.
    ## consider, for example, the person who wants to eat steak for an appetizer, entree, and dessert!
    ##------------------------------------------------------------------------------------------
    PUBLIC$new <- function(candidate_gene_sets,
                           candidate_fan_edges_OR_fan_edges,
                           fan_vertices = NULL,
                           VERBOSE = FALSE)
      {
        if(VERBOSE) print("start: constructing a PrixFixe problem instance.")

        CANDIDATE_GENE_SETS <- candidate_gene_sets
        CANDIDATE_FAN_EDGES <- candidate_fan_edges_OR_fan_edges
        FAN_VERTICES <- fan_vertices


        ##------------------------------------------------------------------------------------------
        ## sanitize FAN_VERTICES.
        ##------------------------------------------------------------------------------------------
        local({
          stopifnot(GPF$ID$NCBI_GENE$isValid(FAN_VERTICES))
          FAN_VERTICES <<- unique(GPF$ID$NCBI_GENE$toFullForm(FAN_VERTICES))
          if(VERBOSE) print(sprintf("%d sanitized FAN vertices", length(FAN_VERTICES)))
        })


        ##------------------------------------------------------------------------------------------
        ## sanitize CANDIDATE_GENE_SETS.
        ##------------------------------------------------------------------------------------------
        local({
          ## should be a list.
          stopifnot(is.list(CANDIDATE_GENE_SETS))

          ## each entry should be a character vector.
          stopifnot(sapply(CANDIDATE_GENE_SETS, is.character))

          ## remove duplicates within each set.
          if(any(sapply(CANDIDATE_GENE_SETS, function(x) any(duplicated(x))))) {
            warning("duplicate genes found within at least one candidate set, uniquifying...", immediate. = TRUE)
            CANDIDATE_GENE_SETS <<- lapply(CANDIDATE_GENE_SETS, unique)
          }

          ## make sure each gene looks like a valid NCBI Gene ID, then force each to the full-form ID.
          stopifnot(sapply(CANDIDATE_GENE_SETS, function(x) all(GPF$ID$NCBI_GENE$isValid(x)))) ## make sure each gene looks like a valid NCBI Gene ID.
          CANDIDATE_GENE_SETS <<- lapply(CANDIDATE_GENE_SETS, function(x) GPF$ID$NCBI_GENE$toFullForm(x))

          ## how many combinations?
          foo <- sapply(CANDIDATE_GENE_SETS, length)
          if(VERBOSE) print(sprintf("candidate gene set sizes = %s", paste0(foo, collapse = ", ")))
          if(VERBOSE) print(sprintf("%G candidate combinations (over %d gene sets).", prod(foo[foo > 0]), length(CANDIDATE_GENE_SETS)))

          ## name the gene sets if not already named.
          if(is.null(names(CANDIDATE_GENE_SETS))) {
            warning("unnamed candidate gene sets, naming automatically as \"set:1\", \"set:2\", ...", immediate. = TRUE)
            names(CANDIDATE_GENE_SETS) <<- sprintf("set:%d", seq_along(CANDIDATE_GENE_SETS))
          }

          ## in case of duplicate set names, uniquify.
          if(any(duplicated(names(CANDIDATE_GENE_SETS)))) {
            warning("duplicated set names found for candidate gene sets, uniquifying...", immediate. = TRUE)
            names(CANDIDATE_GENE_SETS) <<- GPF$UTIL$makeUnique(names(CANDIDATE_GENE_SETS))
          }
        })


        ##------------------------------------------------------------------------------------------
        ## sanitize CANDIDATE_FAN_EDGES, which must be a data frame or matrix with _exactly_ two columns.
        ##
        ## IMPORTANT! also trim to only those genes in CANDIDATE_GENE_SETS.
        ## (in cast a full graph is passed, for example.)
        ##------------------------------------------------------------------------------------------
        local({
          stopifnot(is.data.frame(CANDIDATE_FAN_EDGES) || is.matrix(CANDIDATE_FAN_EDGES))

          stopifnot(ncol(CANDIDATE_FAN_EDGES) >= 2)
          CANDIDATE_FAN_EDGES <<- as.data.frame(CANDIDATE_FAN_EDGES[,1:2], stringsAsFactors = FALSE)

          stopifnot(all(GPF$ID$NCBI_GENE$isValid(c(CANDIDATE_FAN_EDGES[[1]], CANDIDATE_FAN_EDGES[[2]]))))
          CANDIDATE_FAN_EDGES[[1]] <<- GPF$ID$NCBI_GENE$toFullForm(CANDIDATE_FAN_EDGES[[1]])
          CANDIDATE_FAN_EDGES[[2]] <<- GPF$ID$NCBI_GENE$toFullForm(CANDIDATE_FAN_EDGES[[2]])

          candidate_genes <- unique(unlist(candidate_gene_sets))
          lix <- apply(CANDIDATE_FAN_EDGES, 1, function(x) all(x %in% candidate_genes))
          if(any(!lix)) {
            warning("trimming the specified edges to only include those edges incident to pairs of genes in the candidate gene sets.", immediate. = TRUE)
            CANDIDATE_FAN_EDGES <<- CANDIDATE_FAN_EDGES[lix,]
          }
        })


        ##------------------------------------------------------------------------------------------
        ## create a table that will become the main repository for all of the genes within each window/locus.
        ## even for those genes that aren't part of the solution-search procedure, we'd like to report them as candidates to the searcher.
        ##------------------------------------------------------------------------------------------
        CANDIDATE_TABLE <- DataFrame(ldply(1:length(CANDIDATE_GENE_SETS), function(i) {
          data.frame(set_index = rep(i, length(CANDIDATE_GENE_SETS[[i]])),
                     set_name = rep(names(CANDIDATE_GENE_SETS)[i], length(CANDIDATE_GENE_SETS[[i]])),
                     ncbi_gene_id = CANDIDATE_GENE_SETS[[i]],
                     stringsAsFactors = FALSE)
        }))
        if(any(duplicated(CANDIDATE_TABLE$ncbi_gene_id))) {
          warning("vertex sets are not fully disjoint.", immediate. = TRUE)
        }


        ##------------------------------------------------------------------------------------------
        ## is the gene in the (complete) FAN?
        ## what about after trimming the FAN to the candidate edges only?
        ## and then after reducing further to the n-partite graph?
        ##
        ## example: let's say gene A is in the complete FAN.
        ## but A is not connected to any other candidate gene (inside OR outside of A's locus).
        ## then A is in_FAN, but **not** in_CANDIDATE_FAN.
        ##
        ## example: let's say gene B is in the complete FAN.
        ## furthermore, B is **only** connected to C, where B and C happen to be in the same locus.
        ## once making the n-partite graph, this edge disappears, and B is now disconnected.
        ## B is in_FAN,and B is in_CANDIDATE_FAN, but B is **not** in_N_PARTITE_FAN.
        ##
        ## programming note: foo.edges and foo.sets for some reason really need to be plain data frames (i.e. not of class DataFrame), as the latter somehow forces factors into the merged results, which then make checking for equality screwy afterwards.
        ## once the "in_N_PARTITE_FAN" and any other membership expressions are evaluated, the CANDIDATE_FAN_EDGES and N_PARTITE_FAN_EDGES can be DataFrame(d) for prettier interactive use.
        ##------------------------------------------------------------------------------------------
        CANDIDATE_TABLE$in_FAN <- CANDIDATE_TABLE$ncbi_gene_id %in% FAN_VERTICES
        CANDIDATE_TABLE$in_CANDIDATE_FAN <- CANDIDATE_TABLE$ncbi_gene_id %in% c(CANDIDATE_FAN_EDGES[,1], CANDIDATE_FAN_EDGES[,2])
        N_PARTITE_FAN_EDGES <- (function() {
          foo.edges <- as.data.frame(CANDIDATE_FAN_EDGES)
          colnames(foo.edges) <- c("ncbi_gene_id_1", "ncbi_gene_id_2")
          foo.sets <- IRanges::as.data.frame(CANDIDATE_TABLE[,c("set_name", "ncbi_gene_id")])
          foo <- merge(x = foo.edges, y = foo.sets, by.x = "ncbi_gene_id_1", by.y = "ncbi_gene_id")
          foo <- merge(x = foo, y = foo.sets, by.x = "ncbi_gene_id_2", by.y = "ncbi_gene_id")
          foo <- subset(foo, set_name.x != set_name.y)[,c("ncbi_gene_id_1", "ncbi_gene_id_2")]
          return(foo)
        })()
        CANDIDATE_TABLE$in_N_PARTITE_FAN <- CANDIDATE_TABLE$ncbi_gene_id %in% c(N_PARTITE_FAN_EDGES[,1], N_PARTITE_FAN_EDGES[,2])


        ##------------------------------------------------------------------------------------------
        ## now we construct the n-partite igraph object.
        ## we'll actually include all genes in CANDIDATE_TABLE as vertices, just to make our life with igraph easier.
        ## you see, igraph functions that specify vertices not present in the graph will usually error, which is annoying when you want -- for example -- to find the subgraph based on some limited set of vertices.
        ## first restricting the selection set to only those vertices in the main graph adds complexity to the code, so here we'll add all vertices, even if they are completely disconnected in the graph (i.e. with degree == 0).
        ## we'll also add a 'fake' vertex called "foo" that will be used later for vertex-scoring (but again, needs to be present in the graph for when we refer to it).
        ## vertex "foo" is essentially a forced disconnected vertex, in case there aren't any such vertices already.
        ## NOTE: this requires not having any current vertex named "foo" in the graph, so there's a check to make sure.
        ## we could also just make up a unique string and call the vertex that, but "foo" works for now.
        ##
        ## now is the time to also link a few graph properties about each vertex back to the CANDIDATE_TABLE.
        ## this includes degree and the vertex id within the igraph object.
        ## ***** NOTE: these properties refer the n-partite graph. *****
        ##------------------------------------------------------------------------------------------
        G <- graph.data.frame(N_PARTITE_FAN_EDGES, directed = FALSE, vertices = data.frame(c(unique(CANDIDATE_TABLE$ncbi_gene_id), "foo"), stringsAsFactors = FALSE))
        CANDIDATE_TABLE$vertex_id <- match(CANDIDATE_TABLE$ncbi_gene_id, V(G)$name)
        CANDIDATE_TABLE$degree <- degree(G, CANDIDATE_TABLE$vertex_id)


        ##------------------------------------------------------------------------------------------
        ## the 'simplified' vertex sets are those sets with at least one gene having degree > 0.
        ## plus, only those genes are included in the sets (i.e. those with degree == 0 are removed).
        ## the word "sequence" here refers to the idea of a vertex (or edge) sequence in the igraph library.
        ## (i.e. it has nothing to do with biological, or any other kind of, sequences.)
        ##------------------------------------------------------------------------------------------
        foo <- subset(CANDIDATE_TABLE, degree > 0)
        SIMPLIFIED_VERTEX_SEQUENCE_SETS <- IntegerList(tapply(foo$vertex_id, foo$set_name, function(x) x))
        if(VERBOSE) print(sprintf("%G simplified combinations (over %d gene sets).", prod(sapply(SIMPLIFIED_VERTEX_SEQUENCE_SETS, length)), length(SIMPLIFIED_VERTEX_SEQUENCE_SETS)))


        RV <- list()
        RV$CANDIDATE_TABLE <- CANDIDATE_TABLE
        RV$G <- G
        RV$SIMPLIFIED_VERTEX_SEQUENCE_SETS <- SIMPLIFIED_VERTEX_SEQUENCE_SETS
        class(RV) <- "PrixFixe"

        if(VERBOSE) print("done: constructing a PrixFixe problem instance.")

        return(RV)
      }

    return(PUBLIC)
  }

