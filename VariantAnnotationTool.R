rm(list = ls())
library(vcfR)
library(httr)
library(jsonlite)
vcf <- read.vcfR(/data_files/VCFdata.vcf, verbose = FALSE )

fixed = getFIX(vcf) # fix region of vcf

#SECTION 1
# Type of variation: making of 'variant_mutations' vector for all of 6,977 variants
identify.mutations = function(REF,ALT){
  
  if (length(REF) == length(ALT)){
    "SNP" # small nucleotide polymorphism 
  } else if (length(REF) > length(ALT)){
    "DEL" # deletion
  } else if(length(REF) < length(ALT)){
    "INS" # insertion
  } 
  
}

variant_mutations = vector("numeric", nrow(fixed))
for (i in 1:nrow(fixed)){
  REF = unlist(strsplit(fixed[i,4], split=""))
  ALT_vec = unlist(strsplit(fixed[i,5], split=","))
  mutations = vector("character", length(ALT_vec))
  for (z in 1:length(ALT_vec)){
    ALT = unlist(strsplit(ALT_vec[z], split=""))
    mutations[z] = identify.mutations(REF, ALT)
  }
  variant_mutations[i] = paste(mutations, collapse = "|")
}

# SECTION 2
# Depth of sequence coverage at the site of variation
DP = extract.gt(vcf, element = "DP", as.numeric=TRUE)

# SECTION 3
# Number of reads supporting the variant
AO = extract.gt(vcf, element = "AO", as.numeric=TRUE)

#SECTION 4
# Percentage of reads supporting the variant versus those supporting reference
RO = extract.gt(vcf, element = "RO", as.numeric=TRUE)

perc_ALT = round(AO/DP*100) # percent of read supporting variant
perc_REF = round(RO/DP*100) # percent of reads supporting reference

# matrix of perc_ALT/perc_REF for normal and vaf5
M1 = matrix(0, nrow = nrow(vcf), ncol = 2)
for (z in 1:2){
  for (i in 1:nrow(vcf)){
    M1[i,z] = paste0(toString(perc_ALT[i,z]),"/",toString(perc_REF[i,1]))
  }
}

#SECTIONS 5 AND 6
#API
# API helper functions to retrieve variant info in JSON format
API_HELPER = function(base,endpoint){
  call = paste(base, endpoint, sep = "")
  get_var = GET(call)
  get_var_text = content(get_var, "text")
  fromJSON(get_var_text, flatten = TRUE)
}


# website base for the API
base = "http://exac.hms.harvard.edu/rest/variant/variant/"
M2 = matrix("", nrow = nrow(fixed), ncol = 4)
# for-loop for storing ID of variant, its allele frequency, genes, and transcripts
total.time=proc.time()
for (i in 1:nrow(fixed)){
  #Split in case of multiple ALT per variant
  ALT = unlist(strsplit(fixed[i, "ALT"], ","))
  id_vector = vector("character", length(ALT))
  allele_freq_vector = vector("character", length(ALT))
  genes_vector = vector("character", length(ALT))
  transcripts_vector = vector ("character", length(ALT))
  # for-loop for multiple ALTs
  for (z in 1:length(ALT)){
    endpoint = paste0(fixed[i, "CHROM"],"-", fixed[i, "POS"], "-", 
                      fixed[i, "REF"], "-", ALT[z])                         # endpoint for the API call 
    get_var_json = API_HELPER(base, endpoint)                               # API call
    
    id_vector[z] = endpoint                                                 # full ID of variant stored here
    
    if ("allele_freq" %in% names(get_var_json)){                            
      allele_freq_vector[z] = get_var_json$allele_freq                      # allele_freq retrieved
    } else {
      allele_freq_vector[z] = "Not Found"                                   # or allele_freq not found
    }
    
    if ("genes" %in% names(get_var_json)){
      genes_vector[z] = paste(get_var_json$genes, collapse = "|")           # genes retrieved
    } else {
      genes_vector[z] = "Not Found"                                         # or genes not found
    }
    
    if ("transcripts" %in% names(get_var_json)){
      transcripts_vector[z] = paste(get_var_json$transcripts, collapse = "|") # transcripts retrieved
    } else {
      transcripts_vector[z] = "Not Found"                                     # or transcripts not found
    }
  }
  # results organized into a matrix
  M2[i,1] = paste(id_vector, collapse = "||")                                 
  M2[i,2] = paste(allele_freq_vector, collapse = "||")
  M2[i,3] = paste(genes_vector, collapse = "||")
  M2[i,4] = paste(transcripts_vector, collapse = "||")
}
total.time=proc.time()-total.time

# results are organized into a data frame
output_table = data.frame(variant_mutations, DP, AO, M1, M2[,-1])
colnames(output_table) = c("variant_mutations", "normal(ReadDepth)", "vaf5(ReadDepth)",
                           "normal(ALT_Reads)", "vaf5(ALT_Reads)", 
                           "normal(ALT%vsREF%_Reads)", "vaf5(ALT%vsREF%_Reads)", "allele_freq", "genes", "transcripts" )
rownames(output_table) = M2[,1]

head(output_table)