#' @title "functions for working with GTF files"
#' @author Miguel Casanova, mcasanova@medicina.ulisboa.pt
#' @date February 2023
#' @description  The following script, contains a list of functions that can be used to work with GTF files

############################################################################################################################
library(data.table)
library(tidyverse)
library("zeallot")
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
#'  @param removeUnknownAndUncertain parameter to remove repeats that have an unknown origin (for family or class), 
#'  or an uncertain classification (for family or class). This parameter is not used, by default.
#'  @param removeLowComplexity parameter to remove low complexity repeats. The function will remove these by default
#'  @param order parameter used to order the output GTF file by chromosomes and start positions. By default, the function will not order the GTF files
#' OUTPUT:
#'  @return Unmodified GTF file; per TE type GTF file; per individual TE instance GTF file

makeGTF <- function(TE_Annotation,
                    removeUnknownAndUncertain = FALSE,
                    removeLowComplexity = TRUE, 
                    order = FALSE) {
  
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
    toInteger <- c("swScore", "milliDiv", "milliDel", "milliIns", "genoStart", "genoEnd", "genoLeft", "repStart", "repEnd", "repLeft", "id")
    df <- df %>% 
      mutate_at(toInteger, as.integer)
  } else if (source == "U") {
    # Load an UCSC GTF file, downloaded from the table manager section
    df <- fread(file = TE_Annotation, header = T, sep = "\t", data.table = FALSE)
  } else {
    print("You need to input either UCSC or Repeatmasker for this function to work its magic!")
  }
  
  # Several repeats whose annotation is dubious, are labelled with an "?". Moreover, some repeats have an `Unknown` origin.
  # There's two ways of dealing with this: removing all this repeats from the annotation; assuming that the uncertain classification is right
  # and removing the `?` from the family and/or class annotation. 
  # As such, the user has the ability to remove all these uncertain/unknown repeats (not by default), or 
  # remove the `?` characters from both "repClass" and "repFamily"
  if (removeUnknownAndUncertain) {
    repeatMasker_UCSC <- repeatMasker_UCSC %>% 
      filter(!str_detect(repClass, '\\?')) %>%
      filter(!str_detect(repFamily, '\\?')) %>%
      filter(!str_detect(repClass, "Unknown")) %>%
      filter(!str_detect(repFamily, "Unknown")) 
  } else {
    df$repClass <- gsub("\\?", "", as.character(df$repClass))
    df$repFamily <- gsub("\\?", "", as.character(df$repFamily))
  }
  
  # We can now remove the low complexity repeats and repetitive structural RNA.
  # These will be removed by default, but we can choose to keep them.
  if (removeLowComplexity) {
    lowComplexity <- c("Low_complexity", "Simple_repeat", "rRNA", "scRNA", "snRNA", "srpRNA", "tRNA")
    df <- df[ !grepl(paste(lowComplexity, collapse="|"), df$repClass),]
    df <- df[ !grepl(paste(lowComplexity, collapse="|"), df$repFamily),]
  }

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
    df$genoStart <- as.integer(df$genoStart + 1)
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
  
  # # To make sure that the GTF files are properly ordered, we can run the following command
  if (order) {
    GTF <- GTF %>%
      filter(grepl("chr\\d+$", seqname)) %>% # We first filter only canonical chromosomes (chr1-chr22), excluding chrX, chrY and chrM
      arrange(nchar(seqname), seqname, start) %>% # We filter the above rows using nchar (to allow ordering single and double digits correctly), 
      # followed by chromosome name and start position
      bind_rows(GTF %>% 
                  filter(grepl("chr[a-zA-Z]$", seqname)) %>% # To the above, we add the rows correspoding to chrX, chrY and chrM
                  arrange(match(seqname, c("chrX", "chrY", "chrM")), start)) %>% # And arrange them by name and start position
      bind_rows(GTF %>%
                  filter(grepl("chr\\d+_", seqname)) %>% # We then select non-canonical chromosomes (as they have a `_`), excluding sexual chromosomes
                  separate(seqname, c("seqname1","seqname2", "seqname3"), remove = F) %>% # To better sort the chromosome names, these are broken in 3 (separated by `_`)
                  # The first column is the base chromosome, the 2nd is the assembly, the 3rd the patch status
                  arrange(nchar(seqname1), seqname1, seqname3, seqname2, start) %>% # Chromosomes are ordered by chromosome name, followed by patch status (random, alt, fix)
                  # followed by assembly nomemclature and finally the start position
                  select(-seqname1, -seqname2, -seqname3)) %>% # The three columns are removed.
      bind_rows(GTF %>% 
                  filter(grepl("chr[a-zA-Z]_", seqname)) %>%  # Similar to above, but for sexual chromosomes.
                  separate(seqname, c("seqname1","seqname2", "seqname3"), remove = F) %>%
                  arrange(match(seqname1, c("chrX", "chrY", "chrM")), seqname3, seqname2, start) %>%
                  select(-seqname1, -seqname2, -seqname3)) %>%
      bind_rows(GTF %>% 
                  filter(grepl("chrUn_", seqname)) %>% # Similar to above, but for assemblies without known chromosomal origin.
                  separate(seqname, c("seqname1","seqname2"), remove = F) %>%
                  arrange(seqname1, seqname2, start) %>%
                  select(-seqname1, -seqname2))
    } 
  
  # The first GTF file will be used to study types of TEs (using 'gene_name' as the attribute for featuresCount)
  GTF_typeTE <- GTF
  GTF_typeTE$attribute <- gsub("gene_id \"", "gene_name \"TE_", as.character(GTF_typeTE$attribute))
  # The second GTF file will be used to study individual TEs (using 'gene_name' as the attribute for featuresCount)
  GTF_individualTE <- GTF
  GTF_individualTE$attribute <- gsub("transcript_id \"", "gene_name \"TE_", as.character(GTF_individualTE$attribute))
  
  # We now return a list with the three dataframes that can then be loaded individually.
  return(list(original = GTF, TE_Type = GTF_typeTE, TE_Individual = GTF_individualTE))
}

# To run the above function, I use the 'zeallot' library, to create multiple objects from a single return.
# library(zeallot)
# c(UCSC_hg38_rmsk_unmodified, UCSC_hg38_rmsk_TEClass, UCSC_hg38_rmsk_TEIndividual) %<-% makeGTF("pathToTEAnnotation")
