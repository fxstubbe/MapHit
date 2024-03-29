###
# Pipeline 4 - takes the whole_table (generated by the output blast) and retrieves overlapping contigs (end, start & size) + % of overlapping contig
###

### 1) Load packages###
# ----------------------------------------------------------------------------------------------------------- #

suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))

### 2) Functions ###
# ----------------------------------------------------------------------------------------------------------- #

readFastaRef = function(refFile) {
  row = scan(refFile,what=character(0),sep="\n", quiet = T)
  chars = substr(row,1,1)
  base = chars!=">"
  seq = paste(row[base],collapse="")
  return(toupper(unlist(strsplit(seq,""))))
}

readcontigRef = function(refFile) {
  row = refFile
  chars = substr(row,1,1)
  base = chars!=">"
  seq = paste(row[base],collapse="")
  return(toupper(unlist(strsplit(seq,""))))
}


### 3) Paths & Inputs ###
# ----------------------------------------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)
f.root <- getwd()
f.whole_table <- paste( f.root, "/Tables/Whole_table", sep = "" )
f.summary <- paste( f.root, "/Tables/Summary", sep = "" )

#Reference
ref <- readFastaRef(args[1])
ref_ID <- gsub( ".fasta", "", sub('.*\\/', '', args[1] ) )

#Whole table
setwd(f.whole_table)
tables <- list.files(pattern = "\\.csv")
seq_name <- sapply(tables, function(s) gsub(".final._Whole_table.csv","",s))

### 4) Trimming the table ###
# ----------------------------------------------------------------------------------------------------------- #

cat(paste("echo ..... Making a summary Index File .....",sep = " "), sep = "\n")

#sum.file
sum_index.file = paste( f.summary,"/summary_index.txt" , sep = '' )


#------------------Make the cleaned table
for(i in 1:length(tables)){
  
  
  cat(paste("Processing : ", seq_name[i] , sep = ""), sep = "\n" )
  cat(paste(i, "/", length(tables), sep = ""), sep = "\n" )
  
  #setwd(f.alignment)
  
  table <- fread(as.character(tables[i]), select = c(2, 20:22), header = TRUE, sep = ",")
  table <- table %>% distinct()
  colnames(table) <- c('plasmid', 'start', 'end', 'size')
  
  #Add a col containing the plasmid coverage (by unique contigs)
  table$plasmid_coverage <- sapply(table$size, function(x){ (x/length(ref))*100 })
  
  #Add a col with total plasmid coverage
  table$total_plasmid_coverage <- sum(table$plasmid_coverage)
  
  #Add a col informing on matching
  
  for( row in 1 : nrow( table ) ) {
    cat( paste( as.character(seq_name[ i ]), ref_ID, table$start[row], table$end[row],table$size[row],table$plasmid_coverage[row], 
                table$total_plasmid_coverage[row], sep = '\t'), sep = "\n", file = sum_index.file ,append = TRUE )
  }
}

cat(paste("echo ..... DONE Making a summary Index File .....",sep = " "), sep = "\n")
