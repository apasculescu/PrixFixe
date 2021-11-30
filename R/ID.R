## no GPF dependencies.
## initializes a new "ID" object, which acts as a Class with static methods.
initID <- function()
  {
    #require(stringr)
    PUBLIC <- list()


    PUBLIC$DBSNP <- list()
    PUBLIC$DBSNP$pattern <- "^((dbSNP:rs)|(rs))?([0-9]+)$"
    PUBLIC$DBSNP$isValid <- function(dbsnp_ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) dbsnp_ids <- str_trim(dbsnp_ids)
        grepl(PUBLIC$DBSNP$pattern, dbsnp_ids, ignore.case = ignore_case)
      }
    PUBLIC$DBSNP$toFullForm <- function(ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) ids <- str_trim(ids)
        rv <- sprintf("dbSNP:rs%s", sub(PUBLIC$DBSNP$pattern, "\\4", ids, ignore.case = ignore_case))
        rv[!PUBLIC$DBSNP$isValid(ids, ignore_case)] <- NA
        return(rv)
      }
    PUBLIC$DBSNP$toCanonicalForm <- function(ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) ids <- str_trim(ids)
        rv <- sprintf("rs%s", sub(PUBLIC$DBSNP$pattern, "\\4", ids, ignore.case = ignore_case))
        rv[!PUBLIC$DBSNP$isValid(ids, ignore_case)] <- NA
        return(rv)
      }
    PUBLIC$DBSNP$toRsForm <- function(ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) ids <- str_trim(ids)
        rv <- sub(PUBLIC$DBSNP$pattern, "\\4", ids, ignore.case = ignore_case)
        rv[!PUBLIC$DBSNP$isValid(ids, ignore_case)] <- NA
        return(rv)
      }


    PUBLIC$NCBI_GENE <- list()
    PUBLIC$NCBI_GENE$pattern <- "^(NCBI_Gene:)?([0-9]+)$"
    PUBLIC$NCBI_GENE$isValid <- function(gene_ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) gene_ids <- str_trim(gene_ids)
        grepl(PUBLIC$NCBI_GENE$pattern, gene_ids, ignore.case = ignore_case)
      }
    ## "full form" looks like "NCBI_Gene:69"
    PUBLIC$NCBI_GENE$toFullForm <- function(gene_ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) gene_ids <- str_trim(gene_ids)
        rv <- sprintf("NCBI_Gene:%s", sub(PUBLIC$NCBI_GENE$pattern, "\\2", gene_ids, ignore.case = ignore_case))
        rv[!PUBLIC$NCBI_GENE$isValid(gene_ids, ignore_case)] <- NA
        return(rv)
      }
    ## "canonical form" looks like "69"
    PUBLIC$NCBI_GENE$toCanonicalForm <- function(gene_ids, ignore_case = TRUE, str_trim = TRUE)
      {
        if(str_trim) gene_ids <- str_trim(gene_ids)
        rv <- sub(PUBLIC$NCBI_GENE$pattern, "\\2", gene_ids, ignore.case = ignore_case)
        rv[!PUBLIC$NCBI_GENE$isValid(gene_ids, ignore_case)] <- NA
        return(rv)
      }


    return(PUBLIC)
  }
