###
# Pipeline 5 - Output
###


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(tidyverse))
suppressPackageStartupMessages(library(ggplot2))
suppressPackageStartupMessages(library(ggrepel))
suppressPackageStartupMessages(library(RColorBrewer))

### 1) Paths & Inputs ###
# ----------------------------------------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)

f.root <- getwd()
f.sum<- paste(f.root, "/Tables/Summary/summary_index.txt", sep = "")

table <- fread( f.sum , sep = "\t", header=F, data.table = F, fill = T)
colnames(table) <- c("Strain", "Plasmid", "Start", "End", "size", "Plasmid_coverage", "Total_plasmid_coverage")

write.csv(table, "./Tables/Summary/Summary_index.csv", row.names = F)


### 2) iTOL features
# ----------------------------------------------------------------------------------------------------------- #

cat(paste("echo ..... Checking coverage .....",sep = " "), sep = "\n")

boxplot.data <- table[!duplicated(table$Strain), ]

feature = paste("./Tables/Summary/iTOL_full_feature.txt", sep = "")

for( row in 1 : nrow( boxplot.data ) ) { 
  if(boxplot.data$Total_plasmid_coverage[row] > as.numeric( args[1] ) ) {
    cat( paste( as.character( gsub("_", " ",boxplot.data$Strain[row]) ), "1", sep = "\t"), sep = "\n", file = feature , append = TRUE ) 
  } else { 
    cat( paste( as.character( gsub("_", " ",boxplot.data$Strain[row]) ), "0", sep = "\t"), sep = "\n", file = feature , append = TRUE )
    }
} 

#----------------------------------------------------------------------------------------------------------- #



#Table with candidates

Candidates <- table[ which( table$Total_plasmid_coverage >= as.numeric( args[1] ) ) , ]
write.csv(Candidates, "./Tables/Summary/Candidates_index.csv", row.names = F)
Candidates <- Candidates[!duplicated(Candidates$Strain),]
cat(paste( "# candidates : ", nrow(Candidates), sep = ""), sep = "\n")

feature = paste("./Tables/Summary/iTOL_feature.txt", sep = "")

for( row in 1 : nrow( Candidates ) ) { 
  cat( paste( as.character( gsub("_", " ",Candidates$Strain[row]) ), "1", sep = "\t"), sep = "\n", file = feature , append = TRUE ) 
  } 


