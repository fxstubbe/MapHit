##
# Pipeline 3 - Extracting overlapping matches & whole table matches from BLAST output
##

cat(paste("echo ..... Processing blast output .....",sep = " "), sep = "\n")


suppressPackageStartupMessages(library(data.table))
suppressPackageStartupMessages(library(GenomicRanges))
suppressPackageStartupMessages(library(IRanges))
suppressPackageStartupMessages(library(devtools))

### 1) paths and inputs ###
# ----------------------------------------------------------------------------------------------------------- #

args <- commandArgs(trailingOnly = TRUE)
f.root <- getwd()
f.fasta <- args[1]
f.alignment <- paste( f.root, "/Alignment", sep = "" )

#Number of genomes tested

setwd( f.fasta )
fastq <- list.files(pattern = "contigs.fasta")
seq_name <- unique(sapply(fastq, function(s) gsub("contigs.fasta","",s)))
cat(paste("# of genomes =", length(seq_name), sep = " "), sep = "\n")

#Number of BLAST files
setwd( f.alignment )
blast_files <- list.files(pattern = "fasta")
cat(paste("# of blast files =", length(blast_files), sep = " "), sep = "\n")


### 2) Code ###
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


### 3) Code ###
# ----------------------------------------------------------------------------------------------------------- #
setwd(f.root)

#The whole loop goes over each tested genome
for (sr in 1: length(seq_name)) {
  #Progress tag
  cat(paste("Processing : ",gsub(".final.","",seq_name[sr]),sep = "" ),sep = "\n")
  cat(paste("Task : ",sr, "/",length(seq_name),sep = "" ),sep = "\n")
  
  #Selects the alignment file
  NGS = seq_name[sr]
  NGS.1 <- paste(".*",NGS,"contigs.fasta",sep = "")
  my_files <- blast_files[grep(NGS.1, blast_files)]
  
  #Get a proper name
  my_names <- gsub(".final.contigs.fasta","",my_files)
  my_names <- gsub("alignment_pBSRC1.fa.","",my_names)
  
  #Get the name of the plasmid blast onto the strain (it's in the name)
  my_plasmids <- gsub(".final.contigs.fasta","",my_files)
  plasmid_names <- gsub("\\..*","",my_plasmids)
  plasmid_names <- gsub("alignment_","",plasmid_names)
  
  
  #Extract relevant information from blast output file 
  #It's a list of list
  setwd(f.alignment)
  my_list <- lapply(my_files, function(i){fread( i, sep = "\t", header=F, data.table = F, fill = T)})
  setwd(f.root)
  
  #Let's clean the the list of list
  my_clean_list <- list()
  
  
  for(i in 1:length(my_list)){
    #Grep every line
    my_clean_list[[i]] <- my_list[[i]][c(-grep("# ", my_list[[i]]$V1)),]
    
    if(class(my_clean_list[[i]]) == "character"){
      
      if(length(my_clean_list[[i]]) > 0 ){
        
        #Let's make a dataframe called test where we put relevant informations
        test <- sapply(my_clean_list[[i]], function(s) t(as.data.frame(strsplit(s, split = "\t"))))
        test <- as.data.frame(t(test))
        test$V1 <- as.character(test$V1)
        test$V2 <- as.character(test$V2)
        test$V3 <- as.numeric(as.character(test$V3))
        test$V4 <- as.integer(as.character(test$V4))
        test$V5 <- as.integer(as.character(test$V5))
        test$V6 <- as.integer(as.character(test$V6))
        test$V7 <- as.integer(as.character(test$V7))
        test$V8 <- as.integer(as.character(test$V8))
        test$V9 <- as.integer(as.character(test$V9))
        test$V10 <- as.integer(as.character(test$V10))
        test$V11 <- as.numeric(as.character(test$V11))
        test$V12 <- as.numeric(as.character(test$V12))
        test$V13 <- as.integer(as.character(test$V13))
        test$V14 <- as.character(test$V14)
        rownames(test) <- NULL
        
        #Remove each sub-list (one by plasmid) with the corresponding dataframe
        my_clean_list[[i]] <- test
        
      }
    }
  }
  
  #Check if there is something in list
  if (identical(my_clean_list[[1]], character(0) ) ) {next}
  
  names(my_clean_list) <- my_names
  plasmid_match <- my_clean_list[lapply(my_clean_list, length)>1]
  
  #Naming columns
  if(length(plasmid_match) > 0){
    No_plasmid_match <- my_clean_list[lapply(my_clean_list, length)==1]
    column.names <- c("qseqid", "sseqid", "pident", "length_match", 
                      "mismatch", "gapopen","qstart", "qend", 
                      "sstart", "send", "evalue", "bitscore", "sgi", "stitle" )
    plasmid_match <- lapply(plasmid_match, setNames, column.names)
  }
  
  
  
  #--------- select best eValue score for each genome ---------------
  
  plasmid_match_temp1 <- list()
  
  #make a list of every unique contigs
  uniq_contig <- unique(plasmid_match[[1]]$qseqid)
  
  #Every sublist is for a specific contig where good evalue (< 1e-10) match to the genome are kept (sorted by evalue)
  plasmid_match_temp2 <- list()
  
  for (k in 1 : length( uniq_contig ) ) {
    d <- plasmid_match[[1]][which(plasmid_match[[1]]$qseqid == uniq_contig[k]),]
    uniq_plasmid <- unique(d$sseqid)
    plasmid_match_temp3 <- list()
    
    v <- d[which(d$sseqid == uniq_plasmid[1]),]
    w <- v[order(v$evalue),]
    
    plasmid_match_temp3[[1]] <- w[which(w$evalue < 1e-10),]
    
    plasmid_match_temp3 <- do.call(rbind, plasmid_match_temp3)
    plasmid_match_temp2[[k]] <- plasmid_match_temp3
    
  }
  
  #The list of list is fused (by row) into a single item  
  plasmid_match_temp2 <- do.call(rbind, plasmid_match_temp2)
  
  #For each plasmid is created a sublist (in temp_1) where plasmid_temp2 is put
  plasmid_match_temp1[[1]] <- plasmid_match_temp2
  names(plasmid_match_temp1) <- names(plasmid_match)
  
  #Rename into paslmid_match
  plasmid_match <- plasmid_match_temp1
  
  #--------- Add strand column and invert coordinates for "minus" strand -----------
  
  #One loop iteration by plasmid that had matched (in the pBSRC1 case ... only one)
  for (o in 1: length(plasmid_match)) {
    #Addition of a column "minus" full of NA
    plasmid_match[[o]]$strand <- NA
    #For each row (so match of the gicen plasmid)
    for (i in 1:nrow(plasmid_match[[1]])) {
      #Check on which strand it matched and put it in plasmid_match$strand
      if (plasmid_match[[o]]$sstart[i] < plasmid_match[[o]]$send[i]) {
        plasmid_match[[o]]$strand[i] <- "plus"}
      else{ plasmid_match[[o]]$strand[i] <- "minus"}
    }
    
    my_temp <- plasmid_match[[1]]
    #Position inversion if on the minus strand
    for (j in 1:nrow(plasmid_match[[o]])) {
      if(my_temp$strand[j] == "minus"){
        my_temp$sstart[j] <- plasmid_match[[o]]$send[j]
        my_temp$send[j] <- plasmid_match[[o]]$sstart[j]}
    }
    
    #Get modified Plasmid_match (one more column and position inversion)
    plasmid_match[[o]] <- my_temp
  }
  
  #--------- reference size list
  lgth_ref <- list()
  for(o in 1:length(plasmid_match)){
    
    #Get the reference for every blasted plasmid
    plasmid_references <- args[2]
    ref <- readFastaRef(plasmid_references)
    
    #Get the siwe of the plasmid
    lgth_ref[[o]] <- length(ref)  
  }
  
  
  #--------- Add contig size
  
  #Get the contif file related to the tested strain
  my_files <- paste(args[1],"/", fastq[sr], sep = "")
  
  list.contigs <- list()
  
  #get the contigs
  row = scan(my_files,what=character(0),sep="\n",quiet = T)
  chars = substr(row,1,1)
  base = which(chars==">")
  for(i in 1:length(base)){
    chars = substr(row,1,1)
    base = which(chars==">")
    base = append(base, length(chars)+1)
    genome <- row[base[i]:(base[i+1]-1)]
    gen <- chars[base[i]:(base[i+1]-1)]
    seq = gen!=">"
    seq = paste(genome[seq],collapse="")
    list.contigs[[i]] <- seq
    name <- gsub(" .*", "",row[base[i]])
    name <- gsub(">", "",name)
    names(list.contigs[[i]]) <- name
  }
  
  #Create an Empty matrix 
  contig.size <- matrix(NA, length(list.contigs),1)
  
  #Get the length of each contig
  for(i in 1:length(list.contigs)){
    ref <- readcontigRef(list.contigs[[i]])
    contig.size[i,1] <- length(ref)
  }
  
  #Make a dataframe associating the length and the name of the contig
  contig.size <- as.data.frame(contig.size)
  colnames(contig.size) <- "length"
  for(i in 1: length(list.contigs)){
    contig.size$contig_name[i] <- names(list.contigs[[i]])
  }
  
  
  #Add colums with metrics such as the percentage of plasmid covered by the contig
  for(o in 1:length(plasmid_match)){
    for(k in 1:length(plasmid_match[[o]]$qseqid)){
      #Find the postion of the contig (in contig size) corresponding to the contig in plasmid_match
      x <- which(contig.size$contig_name == plasmid_match[[o]]$qseqid[k])[1]
      #Get its length and add it to plasmid_match$contig_size
      plasmid_match[[o]]$contig_size[k] <- as.numeric(contig.size$length[x])
      #Calculate the percentage of plasmid covered by contig
      plasmid_match[[o]]$percentage_plasmid_coverage[k] <- round(as.numeric((plasmid_match[[o]]$length_match[k])/(lgth_ref[[o]])*100))
      #Calculate the percentage of the contig that matches the plasmid
      plasmid_match[[o]]$percentage_contig_coverage[k] <- round(as.numeric((plasmid_match[[o]]$length_match[k])/(contig.size$length[x])*100))
      #ratio of contig / plasmid (total length)
      plasmid_match[[o]]$ratio_contig_plasmid_size[k] <- round(as.numeric((plasmid_match[[o]]$contig_size[k])/(lgth_ref[[o]])*100))
    }
  }
  
  #--------- add columns (for overlaps)
  for(i in 1 : length(plasmid_match)){
    plasmid_match[[i]]$overl_coord_start <- NA
    plasmid_match[[i]]$overl_coord_end <- NA
    plasmid_match[[i]]$overl_size <- NA
  }
  
  
  #--------- calculate coordinates of overlapping reads
  for(o in 1:length(plasmid_match)){
    o = 1
    ir <- IRanges(start = plasmid_match[[o]]$sstart, end = plasmid_match[[o]]$send)
    ol <- findOverlaps(ir, reduce(ir))
    ol <- as.matrix(ol)
    coord <- as.data.frame(reduce(ir))
    
    plasmid_match[[o]]$overlapping_contigs <- ol[,2]
    plasmid_match[[o]]$contigs_number <- paste("overlap_",ol[,2], sep = "")
    
    for (i in 1:nrow(coord)) {
      plasmid_match[[o]]$overl_coord_start[which(plasmid_match[[o]]$overlapping_contigs == i)] <- coord[i,1]
      plasmid_match[[o]]$overl_coord_end[which(plasmid_match[[o]]$overlapping_contigs == i)] <- coord[i,2]
      plasmid_match[[o]]$overl_size[which(plasmid_match[[o]]$overlapping_contigs == i)] <- coord[i,3]
    }
    
    plasmid_match[[o]]$coverage <- round(sum(as.numeric(coord[,3])))
    plasmid_match[[o]]$percentage_ovelap_coverage <- round(sum(as.numeric(coord[,3]))/lgth_ref[[o]]*100)
    plasmid_match[[o]]$plasmid_size <- lgth_ref[[o]]
    plasmid_match[[o]] <- plasmid_match[[o]][order(plasmid_match[[o]]$contigs_number),]
    
  }
  
  #--------- calculate coordinates of overlapping reads for each contig
  
  for(o in 1:length(plasmid_match)){
    #Getting unique contigs by plasmid
    uniq <- unique(plasmid_match[[o]]$qseqid)
    #Create an empty list
    contig_plasmid_match <- list()
    
    #for each plasmid
    for(j in 1:length(uniq)){
      
      #subset plasmid_match for unique contigs
      contig_temp <- plasmid_match[[o]][which(plasmid_match[[o]]$qseqid == uniq[j]),]
      
      #find overlap coordinates
      ir <- IRanges(start = contig_temp$sstart, end = contig_temp$send)
      #Merge overlapping contig together
      ol <- findOverlaps(ir, reduce(ir))
      ol <- as.matrix(ol)
      #Makes a datadrame for the new contigs
      coord <- as.data.frame(reduce(ir))
      write.csv(coord,paste("./Tables/Overlap_match/",NGS,"_overlaps_match.csv", sep = ""), row.names = F)
      
      #Adds a column saying in what overlapping contig is the unique contig put
      contig_temp$overlapping_blast_results_for_each_contig <- ol[,2]
      
      #Adding empty columns
      contig_temp$overl_blast_results_coord_start <- NA
      contig_temp$overl_blast_results_coord_end <- NA
      contig_temp$overlapping_blast_results_size <- NA
      contig_temp$total_blast_coverage <- NA
      contig_temp$ratio_blast_coverage_contig_size <- NA
      contig_temp$ratio_blast_coverage_plasmid_size <- NA
      
      #For each overlapping contig
      for (i in 1:nrow(coord)) {
        #Add overlapping contig info (two unique contigs that jave merged will have the same added value)
        #Overlaping contigs end and start point 
        contig_temp$overl_blast_results_coord_start[which(contig_temp$overlapping_blast_results_for_each_contig == i)] <- coord[i,1]
        contig_temp$overl_blast_results_coord_end[which(contig_temp$overlapping_blast_results_for_each_contig == i)] <- coord[i,2]
        #Overlapping contig size
        contig_temp$overlapping_blast_results_size[which(contig_temp$overlapping_blast_results_for_each_contig == i)] <- coord[i,3]
      }
      
      #Get the total BLAST size
      contig_temp$total_blast_coverage <- round(sum(as.numeric(coord[,3])))
      
      
      #Get ratio blast_coverage over contig size
      contig_temp$ratio_blast_coverage_contig_size <- round(sum(as.numeric(coord[,3]))/(contig_temp$contig_size[j])*100) #weird
      
      
      contig_temp$ratio_blast_coverage_plasmid_size <- round(sum(as.numeric(coord[,3]))/(contig_temp$plasmid_size[1])*100)
      
      contig_plasmid_match[[j]] <- contig_temp
      
    }
    
    plasmid_match[[o]] <- do.call(rbind,contig_plasmid_match)
  }
  
  
  write.csv(plasmid_match,paste("./Tables/Whole_table/",NGS,"_Whole_table.csv", sep = ""), row.names = F)
  
}

cat(paste("echo ..... DONE Processing blast output .....",sep = " "), sep = "\n")
