\name{getLDRegionsFromSNPs.WS}
\title{getLDRegionsFromSNPs.WS}
\alias{getLDRegionsFromSNPs}
\description{Queries, in this version, the Ensembl and NCBI web services to find genomic regions
  and their associated genes based on HapMap LD windows. Requires an active internet/web connection.}
\usage{getLDRegionsFromSNPs.WS(dbsnp_ids,
  r_square_threshold = 0.5,
  upstream_padding = 100,
  downstream_padding = 100,
  populations = c('1000GENOMES:phase_3:CHB','1000GENOMES:phase_3:CEU'))}
\arguments{
  \item{dbsnp_ids}{
    A character vector of valid dbSNP IDs.
    The default value is the `example_dbsnp_ids` vector which is automatically defined in the package.
  }
  \item{r_square_threshold}{
    A numeric scalar between 0.0 and 1.0.
    The default value is 0.5.
  }
  \item{upstream_padding}{
    An integer >= 0 specifying the padding to be applied upstream of
    transcripts (in kilobases).
    The default value is 100. In this version the maximum between upstream_padding
    and downstream_padding will be actually used as a global window for searches.
  }
  \item{downstream_padding}{
    An integer >= 0 specifying the padding to be applied downstream of
    transcripts (in kilobases).
    The default value is 100. In this version the maximum between upstream_padding
    and downstream_padding will be actually used as a global window for searches
  }
  \item{populations}{
    A vector of characters specifying the populations
    The default value is c('1000GENOMES:phase_3:CHB','1000GENOMES:phase_3:CEU')
  }
}
\details{
  It obtains all the dbsnps and their genomic location on the populations and with the restrictions defined in the input arguments.
  From there it extracts the genes in the genomic regions of interest.
  Performs web service calls on a number of public servers, so
  please use judiciously/kindly! Some of the limitations are described here
  https://rest.ensembl.org/documentation/info/ld_id_get and here
  https://eutils.ncbi.nlm.nih.gov
  A dot (.) is echoed on the console after each attempt to retrieve the ld_regions of a dbsnp.
}
\value{
  A list with the following entries:
  \itemize{
    \item{\strong{invalid_dbsnp_ids}}{
      Any specified IDs that don't appear to be dbSNP IDs.
    }
    \item{\strong{unfound_dbsnp_ids}}{
      Valid-looking dbSNP IDs for which no entry can be found.
    }
    \item{\strong{unmapped_dbsnps}}{
      Valid SNPs that don't seem to have been mapped/genotyped as part
      of the HapMap project (i.e. no LD information available).
    }
    \item{\strong{ld_regions}}{
      A data frame called "ld_regions". It has the needed information for PrixFixe to continue.
      Each entry corresponds to a gene that is cytogenetic address name.
      And each locus is itself a list a few components:
      \itemize{
	\item{\emph{region_name}}{
	  :A character vector specifying the genomic regions for each gene.
	  Genomic coordinates correspond to the maximum boundaries of transcripts found
	  within this region, and these transcripts are thus collapsed
	  into individual genes}
	\item{\emph{id}}{
	  :A character vector with gene
	  ids found in the genomic region of interest.}
	\item{\emph{symbol}}{
	  :A character vector with genes names of the gene ids above.}
	\item{\emph{dbSNPs_count}}{
	  :A numeric vector the count of SNPs ids corresponding to the ids above.}
      } % closes the \itemize under \item{ld_regions}
    } % closes \item{ld_regions}
  } % closes \itemize under \value
  The complete structure of the return value is best explored with the function
  str().
}
\examples{
x <- getLDRegionsFromSNPs.WS(example_dbsnp_ids)
str(x)
print(x$ld_regions)
}
