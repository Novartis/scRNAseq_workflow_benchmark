#####################################################
#
# scRNASeq pipeline functions
#
# PART III: Normalization
# _______________________
#
# This script contains a wrapper function to different normalizations, based on the scran R package.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

#########################################
# Normalizations
# Input:
# - sce = an SCESet
# - method = the method you want to use to normalize. Choices are:
#     TC = total count normalization (multiply this by 10^6 to get CPMs)
#     UQ = upperquartile
#     RLE = relative log-expression, as in DESeq2
#     TMM = trimmed mean of M-values, as in edgeR
#     scran (default) = Lun sum factors, implemented in scran package
# Output:
# - the normalized expression values in the norm_exprs slot of the SCESet
#########################################

normalize_counts = function(sce,method = "scran"){
  #calculate size factors according to method
  switch(method,
         "TC" ={sce = normalizeExprs(sce, method="none",return_norm_as_exprs=T)},
         "RLE" = {sce = normalizeExprs(sce, method="RLE",return_norm_as_exprs=T)},
         "TMM" = {sce = normalizeExprs(sce, method="TMM",return_norm_as_exprs=T)},
         "UQ" = {sce = normalizeExprs(sce, method="upperquartile",return_norm_as_exprs=T)}, 
         "scran" = {clusters = quickCluster(sce)
         sizes = seq(20,100,5)
         if(min(table(clusters))>max(sizes)){
           sce = computeSumFactors(sce,clusters = clusters,sizes=sizes)
         } else{
           message("Clustering of cells failed, using global scaling factors")
           sce = computeSumFactors(sce)
           if(any(sizeFactors(sce) < 0)) {
             warning("Negative size factors generated. Most likely, this is due to some cells having very low total feature counts. Consider using more stringent QC cutoffs.")
           }
         }
         sce = scater::normalize(sce, return_norm_as_exprs=T)}
  )
  return(sce)  
}
