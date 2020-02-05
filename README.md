
# Prototype of a Variant Annotation Tool



Information in VCF file was used to create a table annotating each of 
the variants in the file. The output table includes (1) type of genetic variation,
(2) depth of sequence coverage for each variant site, (3) number of reads supporting 
a variant, (4) Percentage of reads supporting the variant vs those supporting the
reference, (5) allele frequency of a variant, (6) genes of a variant location, 
(7) transcripts of a variant location. Items 5, 6, and 7 were acquired from
Broad Institute ExAC Project API. 



The R package vcfR was used to parse the file. For each ALT of a variant
type of variation (mutation), genes of a variant location, and  transcripts of a 
variant location were retrieved for 6,977 variants. 




