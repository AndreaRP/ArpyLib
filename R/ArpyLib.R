setClass("harmony_gene", slots=list(vec1="numeric", vec2="numeric", r1="numeric", r2="numeric", gs.vec="numeric"))


#' capitalize
#'
#' Capitalizes first letter in string
#' @keywords cappitalize
#' @export
#' @examples
#' capitalize("HELLO")
capitalize <- function(x) {
  s <- paste(toupper(substring(x, 1,1)), tolower(substring(x, 2)),
             sep="", collapse=" ")
  return(s)
}

#' mart_data_import
#'
#' Downloads or loads the biomart annotated reference of genes and transcripts (t2g)
#' and the compatec gene to gene symbol reference (t2g.annot)
#' @keywords mart
#' @export
#' @examples
#' annot <- mart.data.import()
#' t2g <- annot$t2g
#' t2g_annot <- annot$t2g_annot

mart_data_import <- function(){
  library("plyr")
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL"
                         , dataset = "mmusculus_gene_ensembl"
                         , host = 'mar2016.archive.ensembl.org'
  )

  # Subset transcript id, gene names and gene ids from biomart
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id"
                                    , "ensembl_gene_id"
                                    , "external_gene_name"
                                    , "mgi_symbol"
                                    , "mgi_id"
                                    , "mgi_description"
  )
  , mart = mart)

  # Rename the mart downloaded frame to please sleuth
  t2g <- dplyr::rename(t2g
                     , target_id = ensembl_transcript_id
                     , ens_gene = ensembl_gene_id
                     , ext_gene = external_gene_name)

  # Create description object
  t2g_annot <- unique(t2g[,c(2:6)])
  # Group descriptions by gene
  t2g_annot <- ddply(t2g_annot
                   , .(ens_gene, ext_gene)
                   , summarize
                   , mgi_description = paste(toString(paste(mgi_symbol, ": ", mgi_description, sep=""))))
  detach("package:plyr", unload=TRUE)
  return(list(t2g=t2g, t2g_annot=t2g_annot))
}


#' mart_data_import_human
#'
#' Downloads or loads the biomart annotated reference of genes and transcripts (t2g)
#' and the compatec gene to gene symbol reference (t2g.annot)
#' @keywords mart
#' @export
#' @examples
#' annot <- mart.data.import()
#' t2g <- annot$t2g
#' t2g.annot <- annot$t2g.annot
mart_data_import_human <- function(){
  library("plyr")
  mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL"
                         , dataset = "hsapiens_gene_ensembl"
                         , host = 'sep2015.archive.ensembl.org'
  )

  # Subset transcript id, gene names and gene ids from biomart
  t2g <- biomaRt::getBM(attributes = c("ensembl_transcript_id"
                             , "ensembl_gene_id"
                             , "external_gene_name"
                             , "description"
  )
  , mart = mart)

  # Rename the mart downloaded frame to please sleuth
  t2g <- dplyr::rename(t2g
                     , target_id = ensembl_transcript_id
                     , ens_gene = ensembl_gene_id
                     , ext_gene = external_gene_name)

  # Create description object
  t2g.annot <- unique(t2g[,c(2:4)])
  # Group descriptions by gene
  t2g.annot <- ddply(t2g.annot
                   , .(ens_gene, ext_gene)
                   , summarize
                   , description = description)
  # write.csv(t2g.annot, "/data3/arubio/references/t2g.annot_GRCh38.82.csv")
  # write.csv(t2g, "/data3/arubio/references/t2g_GRCh38.82.csv")
  detach("package:plyr", unload=TRUE)
  return(list(t2g=t2g, t2g.annot=t2g.annot))
}

#' gg_cols
#'
#' Return a ggplot palette with n colors
#' @param n number of colors to return.
#' @keywords palette
#' @export
#' @examples
#' cols = gg_cols(6)
gg_cols <- function(n){
  gg_colors <- function(n) {
    hues = seq(15, 375, length = n + 1)
    hcl(h = hues, l = 65, c = 100)[1:n]
  }
  cols = gg_colors(n)
  # plot(rep(1,length(cols)),col=cols, pch=19,cex=2)
  return(cols)
}

#' std_err
#'
#' Calculate the standard error of a vector x
#' @param x numeric vector to test.
#' @keywords sd
#' @export
#' @examples
#' sderr = std.err(x)
std_err <- function(x) sd(x)/sqrt(length(x))

#' p2g_import
#'
#' Downloads or loads the biomart annotated reference of protein to gene (p2g)
#' @keywords mart
#' @export
#' @examples
#' p2g <- p2g.import()
p2g_import <- function(){
    # ANNOTATION (optional)
    mart <- biomaRt::useMart(biomart = "ENSEMBL_MART_ENSEMBL"
                           , dataset = "mmusculus_gene_ensembl"
                           , host = 'mar2016.archive.ensembl.org')
    # Subset transcript id, protein names and protein ids from biomart
    p2g <- biomaRt::getBM(attributes = c("uniprot_swissprot"
                                       , "ensembl_gene_id"
                                       , "external_gene_name"
                                       , "mgi_symbol"
                                       # , "mgi_id"
                                       , "mgi_description"
    )
    , mart = mart)
    # Filter out proteins without description
    p2g <- p2g[which(p2g$uniprot_swissprot != ""),]
    # Rename the mart downloaded frame to please sleuth
    p2g <- dplyr::rename(p2g, target_id = uniprot_swissprot,
                         ens_gene = ensembl_gene_id, ext_gene = external_gene_name)

    # Create description object
    # t2g.annot <- unique(t2g[,c(2:6)])
    # Group descriptions by protein
    # p2g <- ddply(p2g, .(target_id, ext_gene), summarize, mgi_description = paste(toString(paste(ext_gene, ": ", mgi_description, sep=""))))
    # p2g <- ddply(p2g, .(target_id), summarize, mgi_description = paste(toString(paste("(",ext_gene, "): ", mgi_description, sep=""))))
    p2g <- dplyr::ddply(p2g, .(target_id), summarize, ext_gene = paste(toString(paste(ext_gene, sep="")))
                      , mgi_description=paste(toString(paste(mgi_description, sep=""))))
    # saveRDS(p2g, "/data3/arubio/references/p2g_mar2016_archive_ensembl.Rds")
  return(p2g)
}

#' NOTin
#'
#' selects all elements of x that are not in y
#' @export
#' @examples
#' c <- c(which(c %!in% z))
'%NOTin%' <- function(x,y)!('%in%'(x,y))


#' at_least_n
#'
#' Determines if at least n values in a vector are higher than 0.
#' @param vector numeric vector to test.
#' @param n number of elements in vector that must be > 0 for the function to return TRUE.
#' @keywords filter
#' @export
#' @examples
#' keeps the rows in exprs that have at least 6 values > 0 (6 samples are > 0)
#' exprs <- exprs[apply(exprs[,1:nrow(s2c)], 1, function(x) at_least_n(x,6)),]
at_least_n <- function(vector, n){
  cumple <- 0
  for(i in vector){
    if(i>0) cumple <- cumple+1
  }
  if(cumple>=n){
    TRUE
  }else{
    FALSE
  }
}

#' at_least_counts_in_samples
#'
#' Determines if at least "s" samples have at least "c" counts in a vector.
#' @param vector numeric vector to test.
#' @param c minimum number of counts there must be in a sample for the function to return TRUE.
#' @param s minimum number of samples that must have "c" number of counts for the function to return TRUE.
#' @keywords filter
#' @export
#' @examples
#' keeps the rows in exprs that have at least 5 counts in 4 samples
#' exprs <- exprs[apply(exprs[,1:nrow(s2c)], 1, function(x) at_least_counts_in_samples(x, 5, 4)),]
at_least_counts_in_samples <- function (vector, c, s) {
  # vector <- tpm[5,s2c$sample]
  # vector <- tpm[4,s2c$sample]
  # counts <- 5
  # samples <- 4
  
  cumple <- 0
  for (i in vector) {
    if (i > c) 
      cumple <- cumple + 1
  }
  if (cumple >= s) {
    TRUE
  }
  else {
    FALSE
  }
}

#' go_terms
#'
#' Functional annotation of a list of genes (GOprofiler)
#' @param gene.list list of genes to be tested.
#' @param universe background list to test against.
#' @keywords annotation
#' @export
#' @examples
#' gene.GO <- go_terms(gene.list=limma.deg, universe=rownames(norm.counts))
go_terms <- function(gene_list, universe){
  go_data <- gprofiler2::gost(query = gene_list
                              , organism = "mmusculus"
                              , ordered_query = FALSE
                              , multi_query = FALSE
                              , significant = TRUE
                              , exclude_iea = FALSE
                              , measure_underrepresentation = FALSE
                              , evcodes = TRUE
                              , user_threshold = 0.05
                              , correction_method = "fdr"
                              , domain_scope = c("annotated", "known", "custom", "custom_annotated")
                              , custom_bg = universe
                              , numeric_ns = ""
                              , sources = NULL
                              , as_short_link = FALSE
  )
  go_data <- go_data$result
  # Reformat table
  if (nrow(go_data)>0){
    go_data$term_id <- gsub(":", ": ", go_data$term_id)
    go_data$term_name <- gsub("; motif:.*", "", go_data$term_name)
    go_data$term_name_fancy <- apply(go_data, 1, function(x) if (nchar(x["term_name"])>60) paste(substr(x["term_name"], 0, 59), "...", sep="") else x["term_name"])
    go_data <- go_data[, c("term_name", "term_id", "source","term_size"
                           , "intersection_size", "p_value"
                           , "query_size", "intersection", "precision", "term_name_fancy")]
    go_data$neglogpval <- -log10(go_data$p_value)
    go_data <- go_data[order(go_data$p_value),]
    # colnames(go_data) <- c("Term Name", "Term ID", "Domain","Term size"
    #                        , "Overlap",  "subgraph.number", "p.value"
    #                        , "query.size", "Genes", "Precision", "Term Name")
  }
  return(go_data)
}


#' data_summary
#'
#' Obtains the sd and mean grouping by column. Will return a df with the mean and sd by group of
#' the chosen variable.
#' from: STHDA barplot tutorial
#' @param data dataframe with data to summarize
#' @param varname column name for the variable which we want to summarize
#' @param groupnames grouping variables
#' @keywords sd, mean, summarize
#' @export
#' @examples
#' results_list <- data_summary(ToothGrowth, varname="len", groupnames=c("supp", "dose"))
#' summarized_data <- data_summary(all_data, varname = "meta2d_AMP", groupnames = "amp")
data_summary <- function(data, varname, groupnames){
  library("plyr")
  summary_func <- function(x, col){
    c(mean = mean(x[[col]], na.rm=TRUE),
      sd = sd(x[[col]], na.rm=TRUE))
  }
  data_sum<-ddply(data, groupnames, .fun=summary_func,
                  varname)
  data_sum <- rename(data_sum, c("mean" = varname))
  return(data_sum)
}


# go_terms <- function(gene.list, universe){
#   # Query
#   go.data <- gProfileR::gprofiler(gene.list
#                                   , organism = "mmusculus"
#                                   , ordered_query = F
#                                   , significant = T
#                                   , exclude_iea = F
#                                   , underrep = F
#                                   , evcodes = F
#                                   , region_query = F
#                                   , max_p_value = 0.05
#                                   , min_set_size = 0
#                                   , max_set_size = 0
#                                   , min_isect_size = 0
#                                   , correction_method = "fdr"
#                                   , hier_filtering = "none"
#                                   , domain_size = "annotated"
#                                   , custom_bg = universe
#                                   , numeric_ns = ""
#                                   , png_fn = F
#                                   , include_graph = T
#                                   , src_filter = c("GO:BP", "GO:CC", "GO:MF", "KEGG"))
#   # Reformat table
#   if (nrow(go.data)>0){
#     go.data$term.id <- gsub(":", ": ", go.data$term.id)
#     go.data$term.name <- gsub("; motif:.*", "", go.data$term.name)
#     go.data$term.name.fancy <- apply(go.data, 1, function(x) if (nchar(x["term.name"])>60) paste(substr(x["term.name"], 0, 59), "...", sep="") else x["term.name"])
#     go.data <- go.data[, c("term.name", "term.id", "domain","term.size"
#                            , "overlap.size",  "subgraph.number", "p.value"
#                            , "query.size", "intersection", "precision", "term.name.fancy")]
#     go.data$neglogpval <- -log10(go.data$p.value)
#     go.data <- go.data[order(go.data$p.value),]
#     # colnames(go.data) <- c("Term Name", "Term ID", "Domain","Term size"
#     #                        , "Overlap",  "subgraph.number", "p.value"
#     #                        , "query.size", "Genes", "Precision", "Term Name")
#   }
#   return(go.data)
# }