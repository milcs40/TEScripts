# TEScripts

This repository contains different R functions to handle transposable elements (TEs) annotations and create GTF files.

## makeGTF.R

This script contains the ```makeGTF``` function. This function was developed to create GTF files from annotations of transposable elements.
It takes annotation files from either UCSC (UCSC Table Download) or RepeatMasker, pre-processes these tables and generates three different types of GTF files: an unmodified GTF file, a GTF file that allows counting TE families and a GTF file that allows counting individual TE instances.

To run the function, I use the `zeallot` library, that allows creating multiple objects from a single return. 
The function can be called using:
`c(rmsk_unmodified_GTF, rmsk_TEClass_GTF, rms_TEIndividual_GTF) %<% makeGTF("pathToTEAnnotation")
