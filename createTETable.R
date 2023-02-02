############################################################################################################################
#' createTETable
#' @author Miguel Casanova
#' 
#' This function was developed to load an annotation file from UCSC or repeatmasker and create a table with information about repeats.
#' It takes annotation files from either UCSC or Repeatmasker, pre-processes the information and generates 
#' a table containing information about individual repeats, including subfamily, family and class, genomic location and size.
#' Three tables are generated: 1 - A table with basic individual repeat information; 2 - A table with information about individual repeat
#' and statistics about it's family and classes; 3 - a summarized table, with condensed information about repeat subfamilies (and stats).
#' 
#' INPUT:
#'  @param TE_Annotation  annotation file from either UCSC (download table tool) or Repeatmasker (automatic TE annotation of given genome build)
#'  @param removeLowComplexity parameter to remove low complexity repeats. The function will remove these by default
#'  @param order parameter used to order the output TE table files by chromosomes and start positions.
#'  By default, the function will not order the GTF files.
#'  @param assemblyReport assembly report file, containing information about the genome build of interest. 
#'  This needs to be downloaded from the NCBI FTP site. By default, downloads `GCA_000001405.28_GRCh38.p13_assembly_report.txt`.  
#' OUTPUT:
#'  @return Table with individual repeat information; extended table with individual repeat information; subfamily table with genomic statistics

createTETable <- function(TE_Annotation, 
                          assemblyReport = "https://ftp.ncbi.nlm.nih.gov/genomes/all/GCA/000/001/405/GCA_000001405.28_GRCh38.p13/GCA_000001405.28_GRCh38.p13_assembly_report.txt",
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
  
  # Several repeats whose annotation is dubious, are labelled with an "?". Let's remove these from both "repClass" and "repFamily"
  df$repClass <- gsub("\\?", "", as.character(df$repClass))
  df$repFamily <- gsub("\\?", "", as.character(df$repFamily))
  
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
  
  # UCSC tables are 0-based position. This is, the start of the sequences start at 0, and not 1. GTFs have to be 1-based. As such, we need to add
  # 1 to all genome start positions.
  # On the other hand, Repeatmasker annotation file, is already 1-based and doesn't need anything done
  if (source == "U") {
    df$genoStart <- as.integer(df$genoStart + 1)
  } 
  
  # With all this done, we can create our repeat table. 
  TETable <- data.frame(chromosome = df$genoName,
                        start = df$genoStart,
                        end = df$genoEnd,
                        size = df$genoEnd- df$genoStart,
                        strand = df$strand,
                        individualRepeat = df$repName_ind,
                        subfamily = df$repName,
                        family = df$repFamily,
                        class = df$repClass)
  
  # # To make sure that the TE table is properly ordered, we can run the following command
  if (order) {
    TETable <- TETable %>%
      filter(grepl("chr\\d+$", chromosome)) %>% # We first filter only canonical chromosomes (chr1-chr22), excluding chrX, chrY and chrM
      arrange(nchar(chromosome), chromosome, start) %>% # We filter the above rows using nchar (to allow ordering single and double digits correctly), 
      # followed by chromosome name and start position
      bind_rows(TETable %>% 
                  filter(grepl("chr[a-zA-Z]$", chromosome)) %>% # To the above, we add the rows correspoding to chrX, chrY and chrM
                  arrange(match(chromosome, c("chrX", "chrY", "chrM")), start)) %>% # And arrange them by name and start position
      bind_rows(TETable %>%
                  filter(grepl("chr\\d+_", chromosome)) %>% # We then select non-canonical chromosomes (as they have a `_`), excluding sexual chromosomes
                  separate(chromosome, c("chromosome1","chromosome2", "chromosome3"), remove = F) %>% # To better sort the chromosome names, these are broken in 3 (separated by `_`)
                  # The first column is the base chromosome, the 2nd is the assembly, the 3rd the patch status
                  arrange(nchar(chromosome1), chromosome1, chromosome3, chromosome2, start) %>% # Chromosomes are ordered by chromosome name, followed by patch status (random, alt, fix)
                  # followed by assembly nomemclature and finally the start position
                  dplyr::select(-chromosome1, -chromosome2, -chromosome3)) %>% # The three columns are removed.
      bind_rows(TETable %>% 
                  filter(grepl("chr[a-zA-Z]_", chromosome)) %>%  # Similar to above, but for sexual chromosomes.
                  separate(chromosome, c("chromosome1","chromosome2", "chromosome3"), remove = F) %>%
                  arrange(match(chromosome1, c("chrX", "chrY", "chrM")), chromosome3, chromosome2, start) %>%
                  dplyr::select(-chromosome1, -chromosome2, -chromosome3)) %>%
      bind_rows(TETable %>% 
                  filter(grepl("chrUn_", chromosome)) %>% # Similar to above, but for assemblies without known chromosomal origin.
                  separate(chromosome, c("chromosome1","chromosome2"), remove = F) %>%
                  arrange(chromosome1, chromosome2, start) %>%
                  dplyr::select(-chromosome1, -chromosome2))
  } 
  
  # We next can create an extra table that gives us a bit more of info about the different repeats. For this, we will add information about
  # subfamilies, families and classes, including the number of elements of each, and the percentage of the genome they occupe.
  # For this, we need to start by loading the assembly report, to calculate the total genome size
  GRCh38p13 <- fread(assemblyReport,
                     skip = "# Sequence-Name",
                     col.names = c("Ensembl","SequenceRole","AssignedMolecule","AssignedMoleculeLocation/Type",
                                   "NCBI", "Relationship", "RefSeq", "Assembly Unit",
                                   "SequenceLenght", "UCSC"))
  
  # As we are working with the primary build of Hg38, we filter to select only primary assembly, plus mitochondrial DNA.
  GRCh38p13Filtered <- GRCh38p13 %>% 
    filter(`Assembly Unit` == "Primary Assembly" | `Assembly Unit` == "non-nuclear")
  # We can then use this to calculate the total size of the hg38 primary genome build.
  genomeSize <- sum(GRCh38p13Filtered$SequenceLenght)
  
  TETable_info <- TETable %>%
    group_by(subfamily) %>%
    mutate(countSubfamily = n(), .after = subfamily) %>% # We will start by counting elements of each subfamily
    mutate(sizeSubfamily = sum(size), .after = countSubfamily) %>% # Followed by estimating their total genomic size
    mutate(sizeSufamilyPerc = round(sizeSubfamily/genomeSize*100, 2), .after = sizeSubfamily)  %>% # And finally, the percentage of the genome they occupy
    #ungroup() %>%
    group_by(family) %>%
    mutate(countFamily = n(), .after = family) %>% # We do the same as above, for families
    mutate(sizeFamily = sum(size), .after = countFamily) %>%
    mutate(sizeFamilyPerc = round(sizeFamily/genomeSize*100, 2), .after = sizeFamily)  %>%
    #ungroup() %>%
    group_by(class) %>%
    mutate(countClass = n(), .after = class) %>% # And finally, for classes of repeats
    mutate(sizeClass = sum(size), .after = countClass) %>%
    mutate(sizeClassPerc = round(sizeClass/genomeSize*100, 2), .after = sizeClass)  %>%
    ungroup()
  
  TESubfamilyStats <- TETable_info %>% 
    dplyr::select(-chromosome, -start, -end, -strand, -individualRepeat, -size) %>% 
    distinct()
  
  # We finally return a list with the two dataframes we just created, and that can be load individually.
  return(list(TE_list = TETable, TE_list_extended = TETable_info, TE_stats = TESubfamilyStats))
}
