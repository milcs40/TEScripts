# TEScripts

This repository contains different R functions to handle transposable elements (TEs) annotations and create GTF files.

## makeGTF.R

This script contains the ```makeGTF``` function. This function was developed to create GTF files from annotations of transposable elements.
It takes annotation files from either UCSC (UCSC Table Download) or RepeatMasker, pre-processes these tables and generates three different types of GTF files: 
- An unmodified GTF file;
- A GTF file that allows counting TE families;
- A GTF file that allows counting individual TE instances.

To run the function, I use the `zeallot` library, that allows creating multiple objects from a single return. 
The function can be called using:
`c(rmsk_unmodified_GTF, rmsk_TEClass_GTF, rms_TEIndividual_GTF) %<% makeGTF("pathToTEAnnotation")`

## createTETable.R

This script contains the ```createTETable``` function. This function was developed to load an annotation file from UCSC or repeatmasker and create a table with information about repeats. It takes annotation files from either UCSC or Repeatmasker, pre-processes the information and generates a table containing information about individual repeats, including subfamily, family and class, genomic location and size. 

Three tables are generated: 
- A table with basic individual repeat information; 
- A table with information about individual repeat and statistics about it's family and classes; 
- A summarized table, with condensed information about repeat subfamilies (and stats).

To run the function, I use the `zeallot` library, that allows creating multiple objects from a single return. 
The function can be called using:
`c(TETable, TETable_Extended, TETable_subfamilyStats) %<-% 
  createTETable("pathToTEAnnotation", 
                removeLowComplexity = F,
                order = T)`

## extract_attributes.R

This script contains the `extract_attributes` function, that I copied from [Biostars](https://www.biostars.org/p/272889/#418833). 
It's a usefull function that allows extracting, from a GTF attributes column, a given attribute.

## chromAnnotConv.R

This script will use a chromosome alias table from [UCSC](http://hgdownload.soe.ucsc.edu/goldenPath/hg38/database/chromAlias.txt.gz).
to convert the chromosome annotation in a GTF file (first column), to a choosen annotation: `UCSC`, `Ensembl`, `Genbank` or `RefSeq`.

## orderGTF_UCSCAnnot.R

This function will take a GTF loaded in memory and reorganize rows according to chromosome and start sites. At the moment, it only orders GTFs with UCSC annotations (e.g. GTFs loaded from UCSC table download, or Repeatmasker).
