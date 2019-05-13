##
# Pipeline 1 - Make a database
##

args <- commandArgs(trailingOnly = TRUE)

Script = "makedb.sh"
cat(paste("echo  " ),sep="\n",file=Script,append=TRUE)
cat(paste("echo ..... Building a Database ....." ),sep="\n",file=Script,append=TRUE)
cat(paste("#!/bin/bash"),sep="\n",file=Script,append=TRUE)
cat( paste( "makeblastdb -in ", as.character(args[1]), 
            " -dbtype nucl" ,
            " -parse_seqids" ,
            " -out ./Database/database_", sep = "" ), sep = "\n", file = Script , append = TRUE)
cat(paste("echo  " ),sep="\n",file=Script,append=TRUE)
cat(paste("echo  " ),sep="\n",file=Script,append=TRUE)
cat(paste("echo " ),sep="\n",file=Script,append=TRUE)
cat(paste("echo ..... DONE ....." ),sep="\n",file=Script,append=TRUE)