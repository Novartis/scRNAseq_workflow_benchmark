#####################################################
#
# scRNASeq pipeline functions
#
# PART VII: Differential expression analysis
# _______________________
#
# This script contains wrapper functions to different differential expression analysis methods.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

##########################################
# Finding marker genes
#########################################

# Inputs for the different functions:
# - sce = the SCESet after QC or after normalization
# - cl_id : the name of the grouping (e.g. "mclust") used in the comparison. Must be
#   a column name of pData(sce)
# - cl_ref: the entry in pData(sce)$cl_id that should be used as the reference,
#   i.e. against which all the rest of the cells are compared
# - alpha = false discovery rate cutoff, default = 0.05
# - fc_cutoff = log2 fold change cutoff, default = 0.5

#______________________________________________
# Wilcoxon test
# Additional input: pseudocount = a constant added to each expression value
# This is only used for calculating fold changes. It wil shrink the fold change
# for very low expressed genes and prevent getting ridiculous values due to division by very
# small numbers. The default is 2, meaning that it will mostly affect genes having
# mean log-transformed expression values below 1.
#______________________________________________

run_wilcoxon_test = function(sce,cl_id,cl_ref,
                             alpha=0.05,fc_cutoff=0.5,pseudocount=NULL){
  
  if(!is.null(pseudocount)){warning("The pseudocount argument is deprecated and will be removed!")}
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  in_clust = pData(sce)[,cl_id] %in% cl_ref
  pvals = apply(norm_exprs(sce),1,
                function(x) wilcox.test(x=x[in_clust],y=x[!in_clust])$p.value)
  fc = apply(norm_exprs(sce),1,
             function(x) mean(x[in_clust])-mean(x[!in_clust]))
  #log2((mean(2^x[in_clust]-1)+pseudocount)/(mean(2^x[!in_clust]-1)+pseudocount))
  
  out = data.table(gene_id = rownames(sce), pval = pvals, log2fc = fc)
  out[,adj_pval:=p.adjust(pval,method="fdr")]
  out[,DE_flag:=(adj_pval<alpha & abs(log2fc)>fc_cutoff)]
  
  return(out)
}
#_________________________
# limma (can be unreliable with zero-inflated data)
# as input for limma-voom, the raw, non-logged counts are provided along with the size factros. 
# The voom transformation then uses tehse to normalize internally.
# as input for limma-trend, log2(CPM+1) are used as recommended by the authors
#____________________
#
run_limma = function(sce,cl_id,cl_ref,
                    alpha=0.05,fc_cutoff=0.5,count_thr=1,pct=50, method = "trend"){
  
  if(!method %in% c("voom","trend")){
    stop("Method has to be either \"voom\" or \"trend\".")
  }
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  res = as.factor(pData(sce)[,cl_id] %in% cl_ref)
  design = model.matrix(~0+res)
  colnames(design) = c("Rest","Cluster")
  rownames(design) = colnames(sce)
  
  # filter out genes not detected at a count of count-thr in at least
  # 50% of cells in at least one cluster. Outlier cells (clusters with only one cell) are ignored.
  clust_sizes = table(pData(sce)[,cl_id])
  clusts = names(clust_sizes[which(clust_sizes>1)])
  
  keep_mat = matrix(rep(NA,dim(sce)[1]*length(clusts)),ncol = length(clusts))
  for(i in seq(length(clusts))){
    keep_mat[,i] = rowSums(counts(sce)[,pData(sce)[,cl_id]==clusts[i]]>=count_thr)>=pct/100*length(which(pData(sce)[,cl_id]==clusts[i]))  
  }
  keep = apply(keep_mat, 1, function(x) any(x))
  
  #convert to DGEList and filter
  dge = convertTo(sce[keep,],"edgeR")
  contrast_matrix = limma::makeContrasts(Cluster-Rest, levels = design)
  
  if(method == "voom"){
    # transforming counts
    voomt = limma::voom(dge,plot=T,design = design)
    
    #do differential expression analysis on voom transformed data
    fit = limma::lmFit(voomt, design)
    fit2 = limma::contrasts.fit(fit,contrast_matrix)
    fit2 = limma::eBayes(fit2)
  } else {
    logCPM = edgeR::cpm(dge, log=TRUE, prior.count=1)
    fit = limma::lmFit(logCPM, design)
    fit2 = limma::contrasts.fit(fit,contrast_matrix)
    fit2 = limma::eBayes(fit2, trend = T)
  }
  
  diff_exp = limma::topTable(fit2,adjust="BH",number = dim(dge)[1])
  out = data.table(gene_id = rownames(diff_exp),diff_exp)
  out[,DE_flag:=as.factor(adj.P.Val < alpha & abs(logFC) > fc_cutoff)]
  
  return(out)
}

#______________________________________________________________________
# MAST
# as input, use either the raw counts as log2(counts+1), stored in the exprs slot
# or the normalized counts on log scale, stored in the norm_exprs slot
#_______________________________________________________________________

run_MAST = function(sce,cl_id,cl_ref,n_cores = 8,nbins=10,min_per_bin=30,
                    alpha=0.05,fc_cutoff=0.5,norm=F,set_thresh=T){
  
  library(MAST)
  
  ignore = is.na(pData(sce)[,cl_id])
  sce = sce[,!ignore] 
  
  options(mc.cores = n_cores)
  if(norm){
    sca = FromMatrix(norm_exprs(sce), pData(sce), fData(sce))} else {
      sca = FromMatrix(log2(counts(sce)+1), pData(sce), fData(sce))
    }
  # adaptive thresholding
  # note how the threshold moves with median expression
  if(set_thresh){
    message("Calculating expression thresholds...\n
            Check the MAST_theresholds plot. If there are no clear bimodal\n
            distributions, the thresholds are likely to be nonsense.\n
            If that is the case, re-run this function setting set_thresh = F")
    thres = thresholdSCRNACountMatrix(assay(sca), nbins = nbins, min_per_bin = min_per_bin)
    if(!any(thres$cutpoint!=0)){message("All cut points are zero. Try using a different
                                        value of nbins and min_per_bin or set set_thresh=FALSE")}
    par(mfrow=c(nbins%%4+1,4))
    plot(thres)
    dev.copy2pdf(file = file.path(plotdir,"MAST_thresholds.pdf"),width=8,height=3*(nbins%%4+1))
    par(mfrow=c(1,1))
    assays(sca) = list(thresh=thres$counts_threshold, counts=assay(sca))
    }
  
  cond=factor(colData(sca)[,cl_id]==cl_ref)
  cond=relevel(cond,"FALSE")
  colData(sca)$condition=cond
  # calculate the cellular detection rate as no. detected features / no. total features
  # and center it 
  colData(sca)$cngeneson = scale(pData(sce)$total_features/dim(sce)[1],scale=F)
  # fit model (note that this will take time for large datasets!!)
  message("Fitting models...")
  zlmCond = zlm(~condition + cngeneson, sca)
  summaryCond = summary(zlmCond, doLRT='conditionTRUE') 
  #extract the results as a data.table
  summaryDt = summaryCond$datatable
  fcHurdle = merge(summaryDt[contrast=='conditionTRUE' & component=='H',.(primerid, `Pr(>Chisq)`)], #hurdle P values
                    summaryDt[contrast=='conditionTRUE' & component=='logFC', .(primerid, coef, ci.hi, ci.lo)], by='primerid') #logFC coefficients
  fcHurdle[,fdr:=p.adjust(`Pr(>Chisq)`, 'fdr')]
  names(fcHurdle)[1:3] = c("gene_id","pval","log2FC")
  fcHurdle[,DE_flag:=as.factor(fdr<alpha & abs(log2FC)>fc_cutoff)]
  return(fcHurdle)
}



