##
# Pipeline 2 - Make BLAST onto the database
##

# Input
args <- commandArgs(trailingOnly = TRUE)
f.root <- getwd()

#Get the contigs
setwd( args[1] )
contigs_fasta <- list.files( pattern = "final.contigs" )
contigs_blast <- gsub( "_contigs.fasta" , "_Blast.txt" , contigs_fasta )
cat( paste("expected blast files : " , length( contigs_blast ) ), sep = "\n")

#Get Plasmids_ID 
signature_ID <- sub('.*\\/', '', args[2] )

#Loop & write the blast command
#Make the Bash file

setwd(f.root)

Script = "blast_contigs_vs_reference.sh"
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
cat(paste("echo  " ),sep="\n",file=Script,append=TRUE)
cat(paste("echo ..... BLAST reference vs contigs ....." ),sep="\n",file=Script,append=TRUE)

for( j in 1 : length( contigs_fasta ) ) {
  
  cat( paste( 'echo Processing : ', contigs_fasta[j], sep = "") , sep = "\n", file = Script , append = TRUE)
  cat( paste( 'echo task : ', j , "/",length( contigs_fasta )," ...", sep = "") , sep = "\n", file = Script , append = TRUE)
  for( i in 1 : 1 ) {
    cat(paste("blastn", 
              " -db ./Database/database_" , 
              " -max_target_seqs 5" , 
              " -outfmt '7 std sgi stitle'" , 
              " -evalue 0.0000000001" , 
              " -query ", args[1] ,"/",contigs_fasta[ j ] , 
              " -out ./Alignment/alignment_" , 
              paste( signature_ID , "." , contigs_blast[ j ] , sep = "" ) , sep = "" ) , sep = "\n", file = Script , append = TRUE )
  }
}

cat(paste("echo ..... BLAST are done ....." ),sep="\n",file=Script,append=TRUE)
cat(paste("echo  " ),sep="\n",file=Script,append=TRUE)
