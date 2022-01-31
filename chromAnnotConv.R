#' chromAnnotConv
#' @author Miguel Casanova
#' 
#' This function will take a GTF file and convert the chromosomes' annotations to a chosen annotation type.
#' The function will automatically load UCSC's `chromAlias.txt` table, and use it for annotation conversion.
#' Annotations can be `UCSC`, `Ensembl`, `Genbank` or `RefSeq`.
#' 
#' INPUT:
#'  @param GTF  GTF table, in a standard format. The script will assume that the chromosome location is in the first column.
#'  @param annotation  string with the annotation type, you want to convert your chromosome column to. Should choose between `UCSC`, `Ensembl`, `Genbank` or `RefSeq`.
#' OUTPUT:
#'  @return Modified GTF table, with all of the chromosomal annotations in column 1, converted to a choosen annotation style.

chromAnnotConv <- function(GTF, annotation, convTable = "http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz") {
  
  # This script will use an annotation conversion table from UCSC (`chromAlias.txt`).
  # Load genome annotation directly from UCSC chromAlias table.
  convTable <- fread(convTable,
                       col.names = c("alias", "UCSC", "source"))
  # The above table is in a long format. Let's cast it into a wide format 
  # Using a tidyverse approach
  convTable <- convTable %>%
    pivot_wider(names_from = source, values_from = alias) %>% #This command will transform the long table, to a wide table format, creating new columns for different sources
    mutate(Ensembl = coalesce(assembly, `assembly,ensembl`)) %>% #The next two lines, will coalesce the Ensembl and Genbank annotations, respectively. 
    mutate(Genbank = coalesce(genbank, `ensembl,genbank`)) %>%
    select(UCSC, refseq, Ensembl, Genbank) %>% #This line will reorder the columns
    rename(UCSC = UCSC, 
           RefSeq = refseq,
           Ensembl = Ensembl,
           NCBI = Genbank)
  
  # # Using a mix of reshape2 and tidy
  # convTable <- dcast(convTable, UCSC ~ source, value.var = "alias")  
  # convTable <- convTable %>%
  #   mutate(Ensembl = coalesce(assembly, `assembly,ensembl`)) %>%
  #   mutate(Genbank = coalesce(genbank, `ensembl,genbank`)) %>%
  #   select(UCSC, refseq, Ensembl, Genbank) %>%
  #   rename(UCSC = UCSC,
  #          RefSeq = refseq,
  #          Ensembl = Ensembl,
  #          NCBI = Genbank)
  
  
  # Next, we will get the index of the column corresponding to the annotation we want to convert to
  if (annotation %in% c("UCSC", "Ensembl", "RefSeq", "NCBI")) {
    annotation <- match(annotation, names(convTable))  
  } else {
    print("Please input one of possible annotation: UCSC, Ensembl, RefSeq and NCBI!")
  }
  
  # Next, we will scan through all the different annotations and match the chromosome column from the GTF file, with
  # any of the four different annotation sources.
  GTFmod <- GTF
  GTFmod$conv <- GTFmod[[1]]
  GTFmod$convEnsembl <- convTable[[annotation]][match(GTFmod[[1]], convTable$Ensembl)]
  GTFmod$convNCBI <- convTable[[annotation]][match(GTFmod[[1]], convTable$NCBI)]
  GTFmod$convRefSeq <- convTable[[annotation]][match(GTFmod[[1]], convTable$RefSeq)]
  GTFmod$convUCSC <- convTable[[annotation]][match(GTFmod[[1]], convTable$UCSC)]
  
  # Assuming that there aren't repeated annotations and only one match will be found,
  # the next line will collect all matches in the `conv` column
  GTFmod <- GTFmod %>%
    mutate(conv = coalesce(convEnsembl, convNCBI, convRefSeq, convUCSC))
  # Using the above chromAlias table, the next command isn't required
  # It will remove strings of `na` or NA entries, and replace them with the name on column 1 of the GTF file 
  # Which would result on no conversion being applied.
  GTFmod$conv <- ifelse(GTFmod$conv == "na" | is.na(GTFmod$conv), GTFmod[[1]], GTFmod$conv)
  # Finally, the column with the conversions will replace the column 1 and the intermediate columns will be deleted.
  GTFmod[[1]] <- GTFmod$conv
  GTFmod <- GTFmod[,-c(10:14)]
  
  return(GTFmod)
}
