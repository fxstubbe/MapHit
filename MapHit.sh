#!/bin/bash


"""

Pipeline : MapHit
Developer : Francois-Xavier Stubbe

The pipeline maps contigs onto a reference (plasmid, phage ...) and returns 
	* overlap_match files : coord of unique contig (fuse overlapping) onto the reference
	* whole_table files : all the informations about the mapping 
	* summary file : infos from whole_table for candidates
	* iTOL_feature for candidates meeting thresehold requirements  

"""
 

### 1) Path(s) and inputs ###
# ----------------------------------------------------------------------------------------------------------- #

Reference="/Users/stubbf02/Fx_Stubbe/ressources/MGE_project/pBSRC1_Mapping/CC8/Brien_MupA/Plasmids/TypeIII_110-517.fasta"
Contigs="/Users/stubbf02/Fx_Stubbe/ressources/genomes/staph/MGE_project/CC8/contigs"
Thresehold=90

#Have to make sure both arguments are given, if not give a warning 
#Have to make sure are existing directory, if not then give a warning


### 2) File Architecture ###
# ----------------------------------------------------------------------------------------------------------- #

mkdir {Alignment,Database} 
mkdir -p Tables/{Summary,Whole_table,Overlap_match}

### 3) Adding Rscripts to the path ###
# ----------------------------------------------------------------------------------------------------------- #


PATH=$PATH:/Users/stubbf02/Fx_Stubbe/ressources/tools/MapHit/


### 4) MGE_Pipeline ###
# ----------------------------------------------------------------------------------------------------------- #

##
# Step 1 : Make a custom BLAST Database
# Arguments : [1] path to reference
##

Rscript MapHit_1.R $Reference
bash makedb.sh

##
# Step 2 : BLAST contigs vs reference
# Arguments : [1] path to contigs [2] path to reference
##

Rscript MapHit_2.R $Contigs $Reference
bash blast_contigs_vs_reference.sh

##
# Step 3 : Extracting overlapping matches & whole table matches from BLAST output
# Arguments : [1] path to contigs [2] path to reference
##

Rscript MapHit_3.R $Contigs $Reference

##
# Step 4 : Create a row by overlapping read (strain, end, start & size) from whole_tables 
#          and makes a summary file
# Arguments : [1] Path to reference
##

Rscript MapHit_4.R $Reference

##
# Step 5 : Create a row by overlapping read (strain, end, start & size) from whole_tables 
#          and makes a summary file
# Arguments : [1] Thresehold
##

Rscript MapHit_5.R $Thresehold

echo ..... ..... ..... ALL DONE, HAVE A GOOD DAY ..... ..... .....