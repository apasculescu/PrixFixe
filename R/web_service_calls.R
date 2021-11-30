getLDRegionsFromSNPs.WS <- function(dbsnp_ids=c("rs498872", "rs2736100", "rs11979158", "rs4295627", "rs6010620", "rs2252586", "rs2157719", "rs4977756", "rs498872", "rs2736100", "rs4295627", "rs2853676", "rs891835", "rs6010620", "rs1412829", "rs6010620", "rs4809324", "rs0", "foo")
                                    , r_square_threshold = 0.5
                                    , upstream_padding = 100
                                    , downstream_padding = 100
                                    , populations=c('1000GENOMES:phase_3:CHB','1000GENOMES:phase_3:CEU')
                                    )
  {

    invalid_dbsnp_ids <- dbsnp_ids[!grepl('^rs[0-9]+$',dbsnp_ids)]
    dbsnp_ids <- unique(setdiff(dbsnp_ids,invalid_dbsnp_ids))

    window <- max(c(downstream_padding,upstream_padding),na.rm=TRUE)

    ##### to get the other related SNPs to a given SNP (using Linkage Disequilibrium LD)
    #####  will use this help: https://rest.ensembl.org/documentation/info/ld_id_get
    ###### we can filter more stringent:
    url_ld <- paste('https://rest.ensembl.org/ld/human/','rs2736100','/1000GENOMES:phase_3:KHV?content-type=application/json',sep='')
    ##### &d_prime=0.5 ( the default value is 0, max value is 1)
    ##### &r2=0.5 ( the default value is 0, max value is 0.85)	#### the r2 limit by Tasan is 0.25 (see http://dalai.mshri.on.ca/prixfixe/GPF/html/index.html)
    ##### & window_size (default value 500  max value 500kb)		#### the mx nmb of bases by Tasan is 50000 (that is 50kb) (see http://dalai.mshri.on.ca/prixfixe/GPF/html/index.html)


    candidate_gene_table <- data.frame()
    url_get_populations <- 'https://rest.ensembl.org/documentation/info/variation_populations'
    message('Searching for SNPs in LD regions')
    for(rs in dbsnp_ids){	### e.g. rs= 'rs2736100'
      for(pop in populations){
        url_ld <- paste('https://rest.ensembl.org/ld/human/',rs,'/',pop,'?d_prime=0.1&r2=',r_square_threshold,'&window_size=',upstream_padding,sep='')
        r_ld <- httr::GET(paste(url_ld))
        ld_info <- httr::content(r_ld, as="text", type="application/json",encoding="UTF-8")
        crt_rs_df <- data.frame(jsonlite::fromJSON(ld_info))
        if('error' %in% colnames(crt_rs_df) ) next;
        if('variation1' %in% colnames(crt_rs_df)){
          candidate_gene_table <- rbind(candidate_gene_table, crt_rs_df)
        }
      }
      #Sys.sleep(0.05);
      cat('.')
    }


    ##### using NCBI to get the SNPs info (if exists and gene name and id and then chromozome map)
    all_rs <- with(candidate_gene_table, unique(c(`variation1`,`variation2` )))

    url_snp_info <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id='
    url_gene_info <- 'https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=gene&id=' #23187' ### this returns a Jason - from which I can extract the maploc "11q23.3"

    ### we will split in max 50 SNP ids
    batch_length <- 50
    rs_list <- split(all_rs, ceiling(seq_along(all_rs)/batch_length))

    all_rs_snp_id   <- vector(mode='character')
    all_rs_gene_name <- vector(mode='character')
    all_rs_gene_id    <- vector(mode='character')
    all_rs_gene_region <- vector(mode='character')

    message('Obtaining Gene info from SNPs')
    k <- 1; # the counter of SNPs
    for(i in 1:length(rs_list)){
      rs <- paste(rs_list[[i]],collapse=',')
      r <- httr::GET(paste(url_snp_info,rs,sep=''))
      snp_info <- httr::content(r, "text",type='text/xml',encoding='UTF-8')
      snp_info_list <- unlist(strsplit(snp_info,"<DocumentSummary ",fixed=TRUE))
      for(j in 1:length(snp_info_list)){
          snp_info_crt <- snp_info_list[j]
          if(snp_info_crt=='') next
          rs_id <- paste('rs',gsub('".+','',gsub('uid="','',snp_info_crt)),sep='')
          gene_name <- gsub('</NAME>.+','',grep('</NAME>',unlist(strsplit(snp_info_crt,'<NAME>')),value=TRUE))
          gene_id	  <- gsub('</GENE_ID>.+','',grep('</GENE_ID>',unlist(strsplit(snp_info_crt,'<GENE_ID>')),value=TRUE))
          if(length(gene_name) != length(gene_id)){
            print(paste( 'Errors in parsing XML result for rs_id from','https://eutils.ncbi.nlm.nih.gov'))
            stop()
          }
          for(n in 1:length(length(gene_name))){
            all_rs_snp_id[k]  <- rs_id[1]
            all_rs_gene_name[k]<- gene_name[n]
            all_rs_gene_id[k]  <- gene_id[n]
            k <- k+1
          }
      }
      Sys.sleep(0.33) #as per NCBI eutils recommendations https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
      cat('.')
    }

    unfound_dbsnp_ids <- setdiff(dbsnp_ids, all_rs_snp_id)

    #### we must also get the chromosome region of the unique gene ids
    ####   some of the rs (SNPs) do not show a gene. for example: https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?db=snp&id=rs2252586
    gene_ids_unique <- unique(all_rs_gene_id[!is.na(all_rs_gene_id)])
    ##### split them similarly in batches of 50 for example
    gene_list <- split(gene_ids_unique, ceiling(seq_along(gene_ids_unique)/batch_length))

    gene_regions_unique <- vector(mode='character')
    gene_id_region_unique <- vector(mode='character')
    message('Obtaining gene regions')
    k=1
    for(i in 1:length(gene_list)){
      gene_ids <- paste(gene_list[[i]],collapse=',')
      g <- httr::GET(paste(url_gene_info,gene_ids,sep=''))
      gene_info <- httr::content(g, "text")
      ##### split by 'Entrezgene ::= '
      gene_info_list <- unlist(strsplit(gene_info,"Entrezgene ::= ",fixed=TRUE))
      for(j in 1:length(gene_info_list)){
        gene_info_crt <- gene_info_list[j]
        if(gene_info_crt=='') next
        geneid_crt <-  gsub(',.+','',gsub('.+ +geneid ','',gene_info_crt))
        gene_region <- gsub('".+','',gsub('.+maploc "','',gene_info_crt))
        gene_regions_unique[k]  <- paste('chr:',gene_region,sep='')
        gene_id_region_unique[k] <- paste('NCBI_Gene:',geneid_crt,sep='')
        k=k+1
      }
      Sys.sleep(0.33) #as per NCBI eutils recommendations https://www.ncbi.nlm.nih.gov/books/NBK25497/#chapter2.Usage_Guidelines_and_Requiremen
      cat('.')
    }

    ###### check unmapped
    unmapped_dbsnps <- all_rs_snp_id[
                            !sapply(all_rs_gene_id
                                   ,function(g) !is.na(g) &
                                     paste('NCBI_Gene:',g,sep='') %in% gene_id_region_unique)
                            ]

    ld_regions <- data.frame(id=gene_id_region_unique
                      , symbol=sapply(gene_id_region_unique
                                    ,function(gid) unique(all_rs_gene_name[paste('NCBI_Gene:',all_rs_gene_id,sep='')==gid]))
                      , region_name=gene_regions_unique
                      , dbSNPs_count = sapply(gene_id_region_unique
                                              ,function(gid) length(unique(all_rs_snp_id[paste('NCBI_Gene:',all_rs_gene_id,sep='')==gid])))
                      , stringsAsFactors=FALSE)

    return(list(
      invalid_dbsnp_ids = invalid_dbsnp_ids,
      unfound_dbsnp_ids = unfound_dbsnp_ids,
      unmapped_dbsnps = unmapped_dbsnps,
#      ld_search_results = candidate_gene_table,
      ld_regions = ld_regions
      ))
  }

getGeneRegionsTable <- function(ld_regions)
  {
    return(as.data.frame(ld_regions))
  }

getCandidateEdges.WS <- function(ncbi_gene_ids, edge_source='MTasan_example')
  {
  if(!(edge_source %in% c('MTasan_example','MTasan_paper','STRING'))) {
    message('ERROR! edge_source must be one of:',paste(c('MTasan_example','MTasan_paper','STRING'),collapse=',') )
    return(NULL)
  }
    #require(httr)
    #require(jsonlite)
    #require(GPF)

    #### for now we return from MTasan_paper
    ### later we could use STRING for example
    if(FALSE){
      main_edges <- read.csv('R/main_FAN_NCBI_ids.csv')
      main_edges_sorted <- main_edges[with(main_edges, NCBI_Gene1>NCBI_Gene2),]
      main_edges_rev <- main_edges[with(main_edges, NCBI_Gene2>NCBI_Gene1),c(2,1)]
      colnames(main_edges_rev) <- colnames(main_edges_sorted)
      main_edges_sorted <- rbind(main_edges_sorted, main_edges_rev)
      with(main_edges_sorted, sum(NCBI_Gene1<NCBI_Gene2))
      with(main_edges_sorted, sum(NCBI_Gene1>NCBI_Gene2))
      main_edges_sorted <- unique(main_edges_sorted) ##### just verifying that

      usethis::use_data(main_edges)

    }
    if(FALSE){
      edges_Example <- RV
      usethis::use_data(edges_Example)
    }

    if(edge_source=='MTasan_example'){
      data(edges_Example)
      RV <- edges_Example
    }
    if(edge_source=='MTasan_paper') {
      data(main_edges)
      ncbi_gene_ids_numeric <- as.numeric(gsub('NCBI_Gene:','',ncbi_gene_ids))
      RV <- main_edges[with(main_edges, `NCBI_Gene1` %in% ncbi_gene_ids_numeric
                          & `NCBI_Gene2` %in% ncbi_gene_ids_numeric),]
    }
    if(edge_source=='STRING') {
      url_string <- 'https://string-db.org/api/tsv/network?required_score=400&species=9606&identifiers='
      ncbi_gene_ids_numeric <- as.numeric(gsub('NCBI_Gene:','',ncbi_gene_ids))
      gene_ids_param <- paste(ncbi_gene_ids_numeric,collapse='%0d')
      r <- httr::GET(paste(url_string,gene_ids_param,sep=''))
      s_text <- httr::content(r, "text",type='text/xml',encoding='UTF-8')
      s_df <- read.delim(text=s_text)
      if('score' %in% colnames(s_df)){
        RV <- s_df[, match(c('preferredName_A','preferredName_B'),colnames(s_df))]
        colnames(RV) <- c('NCBI_Gene1','NCBI_Gene2')
      }else{
        message('Error in retrieving STRING network',s_text)
        return(NULL)
      }
    }

    return(RV)
  }
