#####################################################
#
# scRNASeq pipeline functions
#
# PART I: Setup
# _______________________
#
# This script contains heleper fucntions for reading data, gene annotation and loading required libraries.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

################
# Library Calls
################

# note that these are only what is needed throughout the analysis
# as some of the other packages load a LOT of depndencies, you might have
# to unload them / restart R after using them and before being able to
# load any new packages

library_calls = function(){
  library(Rtsne) 
  library(ggplot2)
  library(data.table)
  library(scater)
  library(scran)
  library(RColorBrewer)
}

################################
# Reading data & Gene annotation
################################

# Helper function to read count matrix from csv. 

read_data = function(infile,experiment_id){
  counts  =  read.delim(infile, header=T, stringsAsFactors=FALSE)
  # matrix count format
  rownames(counts) = counts$Gene.Id
  counts = counts[,-1]
  colnames(counts) =  paste0(experiment_id,"_C",seq(dim(counts)[2]),sep="")
  return(counts)
}


#----------------------------------
# Annotating genes using ensembldb (the advantage is that it is faster than biomaRt, will
# return exactly one entry per gene including non-coding ones, and uses always
# the same version (regardless of bioconductor version), so this is used
# by default)

get_gene_annotations = function(gene_list,v=F,get_descriptions=T,organism = "human"){
  library(ensembldb)
  if(organism == "human"){
    library(EnsDb.Hsapiens.v79)
    edb = EnsDb.Hsapiens.v79
  } else if(organism == "mouse"){
    library(EnsDb.Mmusculus.v75)
    edb = EnsDb.Mmusculus.v75
  } else {stop("Currently, annotation is only available for organisms human or mouse.")}
  
  my_get = function(x,db,v){
    out = tryCatch({get(x,db)[1]},
                   error = function(cond){
                     if(v) {message(cond)}
                     return("")},
                   finally = {})
    return(out)
  }
  
  gene_info = ensembldb::genes(edb,filter = list(GeneIdFilter(gene_list)))
  gene_info_dt = data.table(gene_id = names(gene_info),
                            chr = as.character(seqnames(gene_info)),
                            symbol = make.unique(gene_info$symbol),
                            gene_biotype = gene_info$gene_biotype)
  
  geneName = data.table(gene_id = gene_list)
  geneName = merge(geneName,gene_info_dt,by="gene_id",all=T)
  
  if(get_descriptions & organism == "human"){
    library(org.Hs.eg.db)
    geneName[,eg:=my_get(gene_id,db=org.Hs.egENSEMBL2EG,v),by="gene_id"]
    geneName[,description:=my_get(eg,db=org.Hs.egGENENAME,v),by="eg"]
  } else if(get_descriptions){
    library(org.Mm.eg.db)
    geneName[,eg:=my_get(gene_id,db=org.Mm.egENSEMBL2EG,v),by="gene_id"]
    geneName[,description:=my_get(eg,db=org.Mm.egGENENAME,v),by="eg"]
  }
  return(geneName)
}

#-----------------------------------
# Annotating genes using biomaRt

get_gene_annotations_biomart = function(gene_list){
  library(biomaRt)
  # Load the organism-specific biomart
  ensembl  =  biomaRt::useEnsembl(
    biomart = 'ensembl', 
    dataset = paste0('hsapiens', '_gene_ensembl'),
    version = 83
  )
  #
  geneName  =  biomaRt::getBM(attributes = c('ensembl_gene_id','hgnc_symbol', 
                                            "entrezgene", "description","chromosome_name"), 
                             filters = 'ensembl_gene_id',
                             values = gene_list, 
                             mart = ensembl)
  
  description = lapply(seq(length(geneName$description)),function(i){
    strsplit(geneName$description[i],"[[]")[[1]][1]
  })
  description = (unlist(description))
  geneName$description = description
  colnames(geneName) = c('gene_id','symbol','eg','description','chr')
  geneName = data.table(geneName)
  setkey(geneName,'gene_id')
  geneName = unique(geneName) #remove duplicate entrez gene identifiers
  geneName[,symbol:=make.unique(symbol)]
  #save(geneName,file="data/output/geneName.RData")
  return(geneName)
}
