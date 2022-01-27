#' extract_attributes
#' @author Someone on the internet
#' 
#' This function was copied from the internet. It allows extracting, from a GTF attributes column, a given attribute
#' 
#' INPUT:
#'  @param gtf_attributes  should be the column, in a GTF data frame, that contains the attribute cell (with info on, transcript_id, gene_name, etc...)
#'  @param att_of_interest  string with the attribute of interest, that should be extracted from the attribute cells (e.g. "gene_name")
#' OUTPUT:
#'  @return list of the attributes of interest that are extracted from a attribute's column.

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

