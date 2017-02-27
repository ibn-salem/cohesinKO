#!/bin/bash

#=======================================================================
#
# This script is suppost do document downlaods of public data in the data 
# folder. It should be excuted from within the /data directory. 
#
#=======================================================================

# set some variables here:
BIN=../bin
mkdir -p ${BIN}

#=======================================================================
# General genome assembly based data and tools from UCSC:
#=======================================================================


#-----------------------------------------------------------------------
# UCSC liftover data and tool
#-----------------------------------------------------------------------

# UCSC liftover chains
mkdir -p UCSC
wget -P UCSC http://hgdownload.cse.ucsc.edu/goldenPath/mm9/liftOver/mm9ToMm10.over.chain.gz
gunzip UCSC/*.gz

# download liftOver tool from UCSC:
wget -P ${BIN} http://hgdownload.cse.ucsc.edu/admin/exe/linux.x86_64/liftOver
chmod u+x ${BIN}/liftOver

#=======================================================================
# Hi-C data from Rao et al 2014 Cell
#=======================================================================
mkdir -p Rao2014


# RAO_CELLS="GM12878_primary+replicate HMEC HUVEC HeLa IMR90 K562 KBM7 NHEK CH12-LX"
RAO_CELLS="CH12-LX"

for CELL in ${RAO_CELLS} ; do

    # download
    wget -P Rao2014 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE63nnn/GSE63525/suppl/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz

    # unzip 
    gunzip Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.gz
        
    # re-format TADs into bed file
	if [ "$CELL" == "CH12-LX" ] 
	then		
	    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
			|cut -f 1-3 \
			> Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
	else
	    tail -n +2 Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt \
			|cut -f 1-3 \
			| sed -e 's/^/chr/' \
			> Rao2014/GSE63525_${CELL}_Arrowhead_domainlist.txt.bed 
	fi
	
done

# convert mouse domains from mm9 to mm10 assembly
${BIN}/liftOver \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed  \
	UCSC/mm9ToMm10.over.chain \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed.mm10.bed \
	Rao2014/GSE63525_CH12-LX_Arrowhead_domainlist.txt.bed.mm9.bed_unmapped.bed


#=======================================================================
# Hi-C data form Dixon et al. 2012
#=======================================================================
mkdir -p Dixon2012

#http://132.239.201.216/mouse/hi-c/mESC.domain.tar.gz
# http://132.239.201.216/mouse/hi-c/cortex.domain.tar.gz

#http://chromosome.sdsc.edu/mouse/hi-c/mESC.domain.tar.gz
#http://chromosome.sdsc.edu/mouse/hi-c/cortex.domain.tar.gz

DIXON_MOUSE="mESC cortex"
for C in ${DIXON_MOUSE} ; do

    # download
    wget -P Dixon2012 http://chromosome.sdsc.edu/mouse/hi-c/${C}.domain.tar.gz    
    # extract and rename
    tar xvfz Dixon2012/${C}.domain.tar.gz -C Dixon2012
    cp Dixon2012/${C}/*combined/total.*combined.domain Dixon2012/mouse.${C}.mm9.bed

    # liftover to mm10
	${BIN}/liftOver \
        Dixon2012/mouse.${C}.mm9.bed \
        UCSC/mm9ToMm10.over.chain \
        Dixon2012/mouse.${C}.mm9.bed.mm10.bed \
        Dixon2012/mouse.${C}.mm9.bed_unmapped.bed

done

#=======================================================================
# TADs in mouse and dog from Rudan et al 2015
#=======================================================================
mkdir -p Rudan2015
wget -P Rudan2015 http://www.cell.com/cms/attachment/2026643989/2045457290/mmc2.xlsx

# #=======================================================================
# # Hi-C interactions (within 2Mb)  from Rudan et al 2015
# #=======================================================================
# wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiCseq%5FSummary%2Etxt%2Egz
# 
# wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Fdog%5Fliver%5Fmerged%5F50000%2Etxt%2Egz
# 
# wget -P Rudan2015 ftp://ftp.ncbi.nlm.nih.gov/geo/series/GSE65nnn/GSE65126/suppl/GSE65126%5FHiC%5Fmouse%5Fliver%5Fmerged%5F50000%2Etxt%2Egz
# 
# gunzip Rudan2015/*.gz

#=======================================================================
# Cohesin KO expression data from Dimitris Polychronopoulos <d.polychronopoulos@imperial.ac.uk>
#=======================================================================
# copied file manually to 
# ICL/RNAseqWTRad21KOMacrophages.RData


