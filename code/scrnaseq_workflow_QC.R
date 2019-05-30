#####################################################
#
# scRNASeq pipeline functions
#
# PART II: Quality control
# _______________________
#
# This script contains all functions called from the Quality Control section of the workflow.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

######################
# Quality Control
######################

#_____________
# Gene filters
#___________

# Filtering out genes that are not expressed at a minimum of min_counts in at least n_th cells
gene_filter_by_feature_count = function(counts, n_th , min_counts = 1){
  discard = rowSums(counts>=min_counts) <= n_th
  message(paste('Flagged', length(which(discard)), 'genes.'))
  return(discard)
}

#______________
# Cell filters
#______________

# Filtering out cells with fewer than n_th detected features
cell_filter_by_feature_count = function(counts, n_th){
  discard = colSums(counts>0) <= n_th
  message(paste('Flagged', length(which(discard)), 'cells.'))
  return(discard)
}

# Filtering out cells with fewer than n_th total UMI counts
cell_filter_by_total_UMI = function(counts, n_th){
  discard = colSums(counts) <= n_th
  message(paste('Flagged', length(which(discard)), 'cells.'))
  return(discard)
}

# Filtering out cells with high mitochondrial gene content
calc_mt_content = function(counts, geneName){
  mt = rownames(counts[rownames(counts) %in% rownames(geneName[geneName$chromosome_name=="MT",]),])
  mt_not = setdiff(rownames(counts),rownames(counts[mt,]))
  counts_mt = counts[mt,]
  mt_genes.amount = 100/colSums(counts)*colSums(counts[mt,])
  return(mt_genes.amount)
}

cell_filter_by_mt_content = function(mt_genes.amount, t){
  discard = mt_genes.amount > t
  message(paste('Flagged', length(which(discard)),'cells.'))
  return(discard)
}

#__________________
# Cell cycle phase annotation
#__________________

#using the scran package, which internally calls cyclone (default)

annotate_cell_cycle = function(sce, organism = "human", gene.names = rownames(sce)){
  if(organism == "human"){
    hs.pairs = readRDS(system.file("exdata", "human_cycle_markers.rds", package="scran"))
    assigned = cyclone(sce, pairs=hs.pairs, gene.names = gene.names)} else if (organism == "mouse"){
      mm.pairs = readRDS(system.file("exdata", "mouse_cycle_markers.rds", package="scran"))
      assigned = cyclone(sce, pairs=mm.pairs, gene.names = gene.names)
    } else {stop("Organism has to be human or mouse.")}
  
  return(assigned)
}

# using annotations from cyclebase (provided as cell.cycle.gene.RData)
# note that this is based on mean expression, so it might be a bad 
# idea to use it on un-normalized data

annotate_cell_cycle_custom = function(log2_counts, organism = "human", input_dir = file.path(code_dir, "input_files")){
  
  # Load list of cell-cycle genes from cyclebase
  load(file.path(input_dir,"cell.cycle.gene.RData"))
  head(cell.cycle.gene)
  
  #plot(sort(cell.cycle.gene$periodicity_pvalue))
  cc = cell.cycle.gene[which(cell.cycle.gene$periodicity_pvalue<.05),]
  # 220/361 genes are a subset of the variable gene list of our dataset
  if(organism == "human"){
    cc =  cc[cc$Ensembl.Gene.ID %in% rownames(log2_counts),] 
    
    # selecting genes associetes with cell-cycle phase
    G1.genes = cc[which(cc$peaktime>0 & cc$peaktime<47 ),]$Ensembl.Gene.ID
    S.genes = cc[which(cc$peaktime>47 & cc$peaktime<70 ),]$Ensembl.Gene.ID
    G2.genes = cc[which(cc$peaktime>70 & cc$peaktime<90 ),]$Ensembl.Gene.ID
    M.genes = cc[which(cc$peaktime>90 & cc$peaktime<100 ),]$Ensembl.Gene.ID
    
    G1_S.genes = cc[which(cc$peaktime>40 & cc$peaktime<50 ),]$Ensembl.Gene.ID
    S_G2.genes = cc[which(cc$peaktime>60 & cc$peaktime<70 ),]$Ensembl.Gene.ID
    G2_M.genes = cc[which(cc$peaktime>80 & cc$peaktime<90 ),]$Ensembl.Gene.ID
    
  } else if(organism == "mouse"){
    cc =  cc[cc$Ensembl.Gene.ID.1 %in% rownames(log2_counts),] 
    # selecting genes associetes with cell-cycle phase
    G1.genes = cc[which(cc$peaktime>0 & cc$peaktime<47 ),]$Ensembl.Gene.ID.1
    S.genes = cc[which(cc$peaktime>47 & cc$peaktime<70 ),]$Ensembl.Gene.ID.1
    G2.genes = cc[which(cc$peaktime>70 & cc$peaktime<90 ),]$Ensembl.Gene.ID.1
    M.genes = cc[which(cc$peaktime>90 & cc$peaktime<100 ),]$Ensembl.Gene.ID.1
    
    G1_S.genes = cc[which(cc$peaktime>40 & cc$peaktime<50 ),]$Ensembl.Gene.ID.1
    S_G2.genes = cc[which(cc$peaktime>60 & cc$peaktime<70 ),]$Ensembl.Gene.ID.1
    G2_M.genes = cc[which(cc$peaktime>80 & cc$peaktime<90 ),]$Ensembl.Gene.ID.1
  } else { stop ("Organism has to be either human or mouse")}
  
  
  # Compute the smean of genes associated with each cell-phase
  G1.sum = apply(log2_counts[G1.genes,colnames(log2_counts)],2,mean)
  S.sum = apply(log2_counts[S.genes,colnames(log2_counts)],2,mean)
  G2.sum = apply(log2_counts[G2.genes,colnames(log2_counts)],2,mean)
  M.sum = apply(log2_counts[M.genes,colnames(log2_counts)],2,mean)
  G1_S.sum = apply(log2_counts[G1_S.genes,colnames(log2_counts)],2,mean)
  S_G2.sum = apply(log2_counts[S_G2.genes,colnames(log2_counts)],2,mean)
  G2_M.sum = apply(log2_counts[G2_M.genes,colnames(log2_counts)],2,mean)
  
  # Create a matrix with the sum of genes associated with each cell-cycle pahse per cells
  G1.S.G2.M = rbind(G1.sum[names(sort(G1.sum))],S.sum[names(sort(G1.sum))],
                   G2.sum[names(sort(G1.sum))] ,
                   M.sum[names(sort(G1.sum))],G1_S.sum[names(sort(G1.sum))],
                   S_G2.sum[names(sort(G1.sum))],G2_M.sum[names(sort(G1.sum))]) 
  
  rownames(G1.S.G2.M) = c("G1","S","G2","M","G1_S","S_G2","G2_M")
  
  
  # compute a cell cycle score (divide by totl expression)
  # note that this is biased by expression differences in cc genes,
  # better to use correlations as in Seurat or scran
  cell.phase = do.call(cbind,
                      lapply(seq(dim(G1.S.G2.M)[2]),function(i){
                        res = G1.S.G2.M[,i]/sum(G1.S.G2.M[,i])
                      }))
  colnames(cell.phase) = colnames(G1.S.G2.M)
  cell.phase = t(cell.phase)
  cell.phase = as.data.frame(cell.phase)
  cell.phase = as.data.frame(t(cell.phase),stringsAsFactors = F)
  return(cell.phase)
}

#___________
# QC visualization
#___________

# Plotting QC based on RNA amount detected per cell
plot_RNA_QC = function(input_sce, min_genes, min_UMI){
  par(mfrow=c(1,3))
  hist(log2(input_sce$total_features),xlab="log2[ # detected genes per cell]", main='', cex.axis=1.5,n=100)
  abline(v=min_genes,col=2)
  hist(log2(input_sce$total_counts),xlab="log2 [# of UMIs per cell]", main='', cex.axis=1.5,n=100)
  abline(v=min_UMI,col=2)
  plot(log2(input_sce$total_features),log2(input_sce$total_counts),xlab="log2[ # detected features per cell]",ylab="log2 [total counts per cell]", cex.axis=1.5)
  abline(v=min_genes,col=2)
  abline(h=min_UMI,col=2)
}

# Plotting mitochondrial gene QC
plot_MT_QC = function(sce, t){
  par(mfrow=c(1,1))
  mt_genes.amount = sce$pct_counts_feature_controls_MT 
  #mt gene content per cell
  plot(seq(length(mt_genes.amount)),(mt_genes.amount),pch=10,cex=1.5,col="gray",main=""
       , cex.axis=2,cex.lab=1.5,xlab="cell index",ylab="Ratio of MT-genes[%]")
  points(seq(length(mt_genes.amount))[mt_genes.amount<t],mt_genes.amount[mt_genes.amount<t],pch=10,cex=1.5)
  abline(h=t,col="red",lty=2,cex=5)
  
  #UMI vs no. genes colored by mt gene content
  plotPhenoData(sce, aes_string(x = "log2(total_features)",
                                y = "log2(total_counts)",
                                colour = "pct_counts_feature_controls_MT"))+
    xlab("Total detected features [log2]") + ylab("Total counts [log2]")+
    ggtitle("Total features vs. total counts, colored by MT content")
}

