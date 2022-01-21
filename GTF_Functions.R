#' @title "functions for working with GTF files"
#' @author Miguel Casanova, mcasanova@medicina.ulisboa.pt
#' @date December 2021
#' @description  The following script, contains a list of functions that can be used to work with GTF files

############################################################################################################################
library(data.table)
library(tidyverse)

############################################################################################################################

#' makeGTF
#' @author Miguel Casanova
#' 
#' This function was developed to create GTF files from annotations of transposable elements.
#' It takes annotation files from either UCSC or Repeatmasker, pre-processes the information and generates 
#' three different GTF files from it: a general one (very similar to a normal GTF file), 
#' a per repeat type GTF file and a individual repeat GTF file
#' 
#' INPUT:
#'  @param TE_Annotation  annotation file from either UCSC (download table tool) or Repeatmasker (automatic TE annotation of given genome build)
#' OUTPUT:
#'  @return Unmodified GTF file; per TE type GTF file; per individual TE instance GTF file

makeGTF <- function(TE_Annotation) {
  
  # First of all, we need to determine where was the table downloaded from, as this will lead to different pre-processing steps.
  source <- readline(prompt = "What is the source of your file, UCSC or Repeatmasker?\nPlease input either U (UCSC) or R(Repeatmasker).")
  
  if (source == "R") {
    # If the GTF file is downloaded from Repeatmasker, directly, it will need a little bit different processing.
    df <- fread(file = TE_Annotation, header = FALSE, fill = TRUE)
    # We can remove the first three lines of the table, that have the old header rows and an empty row.
    df <- df[-c(1:3),]
    
    # Let's set the header names:
    header_final <- c("swScore", "milliDiv", "milliDel", "milliIns", "genoName", "genoStart", "genoEnd", "genoLeft", "strand", "repName", "repClass/repFamily", "repStart", "repEnd", "repLeft", "id")
    names(df) <- header_final
    
    # We can now split the "repClass/repFamily, into two. This will give us information about repeat class and repeat family, in two separated columns.
    df2 <- data.frame(do.call('rbind', strsplit(as.character(df$`repClass/repFamily`), '/', fixed=TRUE)))
    names(df2) <- c("repClass", "repFamily")
    
    df <- df %>%
      add_column("repClass" = df2$repClass,
                 "repFamily" = df2$repFamily,
                 .after = "repClass/repFamily") %>%
      select(-11)
    
    # Repeatmasker annotation files have, on the "strand" column, "C" instead of "-". Let's change all "C"s to "-", as this is required for featureCounts
    df$strand <- gsub("C", "-", as.character(df$strand))
    
    # We next need to correct a few columns and transform them into numeric.
    # Certain columns in the repeatmasker file, have "()" for negative positions on the "genoLeft", "repStart" and "repLeft". Let's correct this.
    df$genoLeft <- gsub("\\(", "-", as.character(df$genoLeft))
    df$genoLeft <- gsub("\\)", "", as.character(df$genoLeft))
    df$repStart <- gsub("-", "", as.character(df$repStart))
    df$repStart <- gsub("\\(", "-", as.character(df$repStart))
    df$repStart <- gsub("\\)", "", as.character(df$repStart))
    df$repLeft <- gsub("\\(", "-", as.character(df$repLeft))
    df$repLeft <- gsub("\\)", "", as.character(df$repLeft))
    
    # We can now define the columns we want to convert to numeric
    toNumeric <- c("swScore", "milliDiv", "milliDel", "milliIns", "genoStart", "genoEnd", "genoLeft", "repStart", "repEnd", "repLeft", "id")
    df <- df %>% 
      mutate_at(toNumeric, as.numeric)
  } else if (source == "U") {
    # Load an UCSC GTF file, downloaded from the table manager section
    df <- fread(file = TE_Annotation, header = T, sep = "\t", data.table = FALSE)
  } else {
    print("You need to input either UCSC or Repeatmasker for this function to work its magic!")
  }
  
  # Several repeats whose annotation is dubious, are labelled with an "?". Let's remove these from both "repClass" and "repFamily"
  df$repClass <- gsub("\\?", "", as.character(df$repClass))
  
  df$repFamily <- gsub("\\?", "", as.character(df$repFamily))
  
  # We can now remove the low complexity repeats and repetitive structural RNA.
  lowComplexity <- c("Low_complexity", "Simple_repeat", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA")
  df <- df[ !grepl(paste(lowComplexity, collapse="|"), df$repClass),]
  df <- df[ !grepl(paste(lowComplexity, collapse="|"), df$repFamily),]
  
  # Once the table has been cleaned, we will proceed with counting repeats to have a way of tracing individual ones.
  # This will count each repeat, in a sequential manner, adding the number of the current iterations of the repeat to the row.
  counter <- with(df, ave(as.character(repName), repName, FUN = seq_along))
  
  # Now, we create a new column with the count values.
  df <- df %>%
    add_column("repName_ind" = counter,
               .after = "repName")
  # We are just left with renaming the count column, by adding the repeat name, followed "_dup" and the count value.
  # This will effectively produce a string with individual repeat entries. E.g. "LTR7_dup20"
  df$repName_ind <- paste0(df$repName, "_dup",  df$repName_ind)
  
  # Finally, we need to create an attribute column. This will contain the following attributes:
  str <- 'gene_id "%s"; transcript_id "%s"; family_id "%s"; class_id "%s"'
  # Let's create the attribute column, for all rows:
  attribCol <- sprintf(str, df$repName, df$repName_ind, df$repFamily, df$repClass)
  # Now we only have to add it to the dataframe
  df <- df %>%
    add_column("attributes" = attribCol,
               .after = "id")
  
  # UCSC tables are 0-based position. This is, the start of the sequences start at 0, and not 1. GTFs have to be 1-based. As such, we need to add
  # 1 to all genome start positions.
  # On the other hand, Repeatmasker annotation file, is already 1-based and doesn't need anything done
  if (source == "U") {
    df$genoStart <- df$genoStart + 1
  } 
  
  # With all this done, we can create the GTF files. 
  GTF <- data.frame(seqname = df$genoName,
                    source = "hg38_rmsk", 
                    feature = "exon",
                    start = df$genoStart,
                    end = df$genoEnd,
                    score = df$swScore,
                    strand = df$strand,
                    frame = ".",
                    attribute = df$attributes)
  # The first GTF file will be used to study types of TEs (using 'gene_name' as the attribute for featuresCount)
  GTF_typeTE <- GTF
  GTF_typeTE$attribute <- gsub("gene_id \"", "gene_name \"TE_", as.character(GTF_typeTE$attribute))
  # The second GTF file will be used to study individual TEs (using 'gene_name' as the attribute for featuresCount)
  GTF_individualTE <- GTF
  GTF_individualTE$attribute <- gsub("transcript_id \"", "gene_name \"TE_", as.character(GTF_individualTE$attribute))
  
  # We now return a list with the three dataframes that can then be loaded individually.
  return(list(original = GTF, TE_Type = GTF_typeTE, TE_Individual = GTF_individualTE))
}



# Function to extract a name, from an attribute, from the attribute column of a GTF file.

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- strsplit(gtf_attributes, "; ")
  att <- gsub("\"","",unlist(att))
  if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
    return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
  }else{
    return(NA)}
}

extract_attributes(hg38GTF_TEType$attribute, "transcript_id")[1]

####
extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

############################################################################################################################

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


# Function to extract a name, from an attribute, from the attribute column of a GTF file.

# This function takes longer to run. Use the next one!
# extract_attributes1 <- function(gtf_attributes, att_of_interest){
#   att <- strsplit(gtf_attributes, "; ")
#   att <- gsub("\"","",unlist(att))
#   if(!is.null(unlist(strsplit(att[grep(att_of_interest, att)], " ")))){
#     return( unlist(strsplit(att[grep(att_of_interest, att)], " "))[2])
#   }else{
#     return(NA)}
# }

extract_attributes <- function(gtf_attributes, att_of_interest){
  att <- unlist(strsplit(gtf_attributes, " "))
  if(att_of_interest %in% att){
    return(gsub("\"|;","", att[which(att %in% att_of_interest)+1]))
  }else{
    return(NA)}
}

