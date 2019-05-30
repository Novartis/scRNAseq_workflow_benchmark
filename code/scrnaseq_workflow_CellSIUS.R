#####################################################
#
# scRNASeq pipeline functions
#
# PART VI: Cell subtype identification using CellSIUS
# _______________________
#
# This implements the CellSIUS method. Note that this can also be used with a custom workflow, as long as you provide
# your data as an SCESet that contains all relevant fields in pData and fData.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################

#-----------------------------------------------
#                 LICENSE
#-----------------------------------------------
#Copyright 2018 Novartis Institutes for BioMedical Research Inc.
# Licensed under the Apache License, Version 2.0 (the "License");
# you may not use this file except in compliance with the License.
# You may obtain a copy of the License at
# http://www.apache.org/licenses/LICENSE-2.0
# Unless required by applicable law or agreed to in writing, software
# distributed under the License is distributed on an "AS IS" BASIS,
# WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
# See the License for the specific language governing permissions and
# limitations under the License.
# limitations under the License.



#########################################
# Identifying rare cell types
#########################################



cellsius_main = function(sce,group_id,min_n_cells=10,verbose = T, min_fc = 2,
                                     organism = "human", corr_cutoff = NULL, iter=0, max_perc_cells = 50,
                                     fc_between_cutoff = 1, mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(Ckmeans.1d.dp)
  
  expr_dt = data.table(gene_id = rownames(sce),norm_exprs(sce))
  expr_dt_melt = melt(expr_dt,id.vars="gene_id",val="expr",var="cell_idx")
  expr_dt_melt = merge(expr_dt_melt,
                       data.table(cell_idx=colnames(sce),main_cluster=as.character(pData(sce)[,group_id])),
                       by="cell_idx")
  
  #Identify genes with significant bimodal distribution
  
  expr_dt_melt[,c("N_cells","within_p","pos0","pos1","Dpos"):=cellsius_find_bimodal_genes(expr,min_n_cells = min_n_cells, max_perc_cells = max_perc_cells),by=c('gene_id','main_cluster')]
  expr_dt_melt[,sig := within_p<100 & Dpos > min_fc] 
  expr_dt_melt[sig==T, within_adj_p:=p.adjust(within_p),by=c('cell_idx')] #correct for multiple testing, only consider genes where test has actually been run
  expr_dt_melt[,sig:=within_adj_p<0.1] 
  expr_dt_melt = expr_dt_melt[gene_id %in% expr_dt_melt[!is.na(sig) & sig==T]$gene_id] 
  
  # If no bimodal gene were found, exit and return NA
  if(dim(expr_dt_melt)[1] == 0){
    print("No genes with bimodal distribution found, returning NA.")
    return(NA)
  }
  # Check whether these genes are specific to the subcluster
  
  for(clust in unique(expr_dt_melt$main_cluster)){
    expr_dt_melt = expr_dt_melt[,paste0(clust,"_",c("p_between","fc")):=cellsius_test_cluster_specificity(
      expr,main_cluster,clust, fc_between_cutoff = fc_between_cutoff),by="gene_id"]
    
    expr_dt_melt[main_cluster==clust,keep:=(expr_dt_melt[main_cluster==clust][[paste0(clust,"_p_between")]] < 0.1)]
  }
  
  expr_dt_melt = expr_dt_melt[keep==TRUE & !is.na(sig)]
  
  # If there are still non-specific genes, discard them (this can happen for
  # very high expressed genes like mitochondrial genes)
  expr_dt_melt[,n_clust_per_gene:=length(unique(main_cluster)),by='gene_id']
  expr_dt_melt = expr_dt_melt[n_clust_per_gene==1]
  expr_dt_melt[,n_clust_per_gene:=NULL]
  
  # Identify correlated gene sets with MCL
  expr_dt_melt = expr_dt_melt[,gene_cluster:=0]
  expr_dt_melt = cellsius_find_gene_sets(expr_dt_melt, corr_cutoff = corr_cutoff, mcl_path=mcl_path)
  
  # discard gene sets that only contain one gene (those are assigned to cluster 0)
  expr_dt_melt = expr_dt_melt[gene_cluster !=0 ]
  
  if(dim(expr_dt_melt)[1] == 0){
    print("No subclusters found, returning NA.")
    return(NA)
  }
  
  # Extract cell subclusters
  expr_dt_melt[,sub_cluster:=main_cluster]
  expr_dt_melt[,mean_expr := mean(expr), by = c('main_cluster','gene_cluster','cell_idx')]
  expr_dt_melt[,sub_cluster:=cellsius_sub_cluster(mean_expr,sub_cluster,gene_cluster, iter=iter),by=c('main_cluster','gene_cluster')]
  
  # Check how many cells belong to the subgroup relative to the total cluster size.
  # If a sub cluster contains more than max_perc_cells cells, discard it.
  clust_list = expr_dt_melt[,list(sub = length(unique(cell_idx))) ,by=c('sub_cluster','main_cluster')]
  clust_list[,tot := sum(sub)/(length(sub_cluster)/2), by= 'main_cluster']
  clust_list = clust_list[grep('_1$',sub_cluster)]
  clust_list[,perc:=sub/tot*100]
  discard_sub_clust = clust_list[perc > max_perc_cells]$sub_cluster
  discard_sub_clust = append(discard_sub_clust,gsub('_1$','_0',discard_sub_clust))
  
  expr_dt_melt = expr_dt_melt[!sub_cluster%in%discard_sub_clust]
  
  # If verbose is TRUE, print a summary of the results
  if(verbose){
    # annotate genes (only if verbose)
    gene_info = get_gene_annotations(unique(expr_dt_melt$gene_id),get_descriptions = T,
                                     organism = organism)
    expr_dt_melt = merge(expr_dt_melt,gene_info, by = 'gene_id')
    
    cellsius_print_summary(expr_dt_melt)
  }
  return(expr_dt_melt)
}

##################################################
# STEP 1: Identify genes with bimodal distribution
##################################################

cellsius_find_bimodal_genes = function(expr, min_n_cells, max_perc_cells){
  
  #skip genes with 0 expression
  if(sum(expr)==0){
    return(list(-1,100,-1,-1,-1))
  }
  # run k-means
  k1d = Ckmeans.1d.dp(expr,k=2)
  # check if a cluster with more than n cells exists
  indx = which(k1d$size>min_n_cells) 
  if(length(indx)>1 ){
    
    # do statistic only if in pos2 cells are less than max_perc_cells% of the total cells in the cluster
    if(k1d$size[2] < round(length(expr)*max_perc_cells/100)){ 
      
      t1=tryCatch({t.test(expr[which(k1d$cluster==2)],y=expr[which(k1d$cluster==1)])},
                   error = function(cond){return(0)},
                   finally={}
      )
      
      if(!is.numeric(t1)){
        
        p1=t1$p.value
        N0=k1d$size[1] # number of cells where the gene is downregulated
        N1=k1d$size[2] # number of cells  where the gene is upregulated
        pos0=k1d$centers[1] 
        pos1=k1d$centers[2]
        Dpos=pos1-pos0
        return(list(N1,p1,pos0,pos1,Dpos))
      } #else {print(paste("ttest failed, dpos = ",pos1-pos0))} # for testing
    }
  }
  # if no cluster was found, return a list of dummy values
  return(list(-1,100,-1,-1,-1))
}

##################################################
# STEP 2: Check whether these genes are specific to one cell subgroup
###################################################

cellsius_test_cluster_specificity = function(exprs, cluster, current_cluster, fc_between_cutoff){
  
  in_clust = which(cluster == current_cluster)
  k1d = Ckmeans.1d.dp(exprs[in_clust],k=2)
  in_subclust = in_clust[which(k1d$cluster==2)]
  
  mean_in = mean(exprs[in_subclust])
  mean_out = mean(exprs[-in_subclust])
  mean_out_nozero = mean(exprs[-in_subclust][exprs[-in_subclust]>0])
  
  # If there are subclusters, but all cells outside the subcluster express 0,
  # set mean_out_nozero to 0
  if(length(in_subclust>0) && !any(exprs[-in_subclust]>0)){mean_out_nozero=0}
  
  fc = mean_in - mean_out
  
  ts = tryCatch({t.test(exprs[in_subclust],exprs[-in_clust])},
                error = function(cond){ return(0)})
  
  if(!is.numeric(ts)){pv = ts$p.value} else {
    #print(paste("ttest failed, fc = ",mean_in-mean_out_nozero)) #for testing only
    pv=999}
  
  if(!is.nan(mean_out_nozero) && (mean_in-mean_out_nozero < fc_between_cutoff)) pv = 999
  return(list(pv,fc))
}

#####################################################
# STEP 3: MCL clustering to find correlated gene sets
#####################################################

cellsius_find_gene_sets = function(expr_dt_melt, corr_cutoff = NULL, min_corr = 0.35, max_corr = 0.5,
                          mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(igraph)
  
  for(clust in unique(expr_dt_melt$main_cluster)){
    
    if(length(unique(expr_dt_melt[main_cluster == clust]$gene_id))==1) { next }
    
    mat = dcast.data.table(expr_dt_melt[main_cluster==clust], gene_id ~ cell_idx,
                           value.var = 'expr')
    mat = mat[rowSums(mat[,-1,with=F])!=0,]
    corr.mat = cor(t(mat[,-1,with=F]))
    dimnames(corr.mat) = list(mat$gene_id,mat$gene_id)
    
    if(is.null(corr_cutoff)){
      corr_cutoff = max(quantile(corr.mat[corr.mat!=1],0.95),min_corr)
      corr_cutoff = min(corr_cutoff, max_corr)}
    adj.corr = corr.mat
    adj.corr[adj.corr<corr_cutoff] = 0
    adj.corr[adj.corr>=corr_cutoff] = 1
    diag(adj.corr) = 0 # no self-loop for MCL
    
    graphs = get.data.frame( graph_from_adjacency_matrix(adj.corr), what = "edges") # gene connection for graphs
    
    # if a graph has no edges (i.e. all genes are uncorrelated), 
    # assign all genes to cluster "0" and go to next iteration
    if(dim(graphs)[1]==0){
      expr_dt_melt = expr_dt_melt[main_cluster == clust, gene_cluster := 0]
      next
    }
    
    graphs = data.frame(graphs,CORR=sapply(seq(dim(graphs)[1]), function(i) corr.mat[graphs$from[i],graphs$to[i]] -corr_cutoff))
    write.table(graphs, file = "tmp.mcl.inp",row.names=F,col.names=F,sep = " ")
    system(paste0(mcl_path, " tmp.mcl.inp --abc -o tmp.mcl.out"))
    x = scan("tmp.mcl.out", what="", sep="\n")
    y = strsplit(x, "[[:space:]]+")
    y = lapply(seq(length(y)), function(i){
      tmp = sapply(seq(length(y[[i]])),function(j){
        gsub('\"','',y[[i]][j])
      })
    })
    
    for(i in seq(length(y))){
      if(length(y[[i]] > 1)){
        expr_dt_melt = expr_dt_melt[main_cluster==clust & gene_id %in% y[[i]],gene_cluster:=i]
      }
    }
  }
  
  return(expr_dt_melt)
}  

############################################
# Step 4: Assign cells to subgroups
############################################

cellsius_sub_cluster = function(mean_expr,sub_cluster,gene_cluster, iter = 0){
  
  k1d = Ckmeans.1d.dp(mean_expr,k=2)$cluster
  cells_sub = (k1d==2)
  
  if(iter == 0){return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))}
  
  # if iter is set higher than 0, a second step of kmeans clustering. 
  # This will remove the lowest peak and can sometimes help to get a more
  # accurate classification.
  
  k1d = Ckmeans.1d.dp(mean_expr[cells_sub],k=2)$cluster
  
  if (max(k1d)>1) {
    cells_sub[cells_sub] = (k1d==2)
    return(paste0(sub_cluster,"_",gene_cluster,"_",as.numeric(cells_sub)))
  }
  return(paste0(sub_cluster,"_",gene_cluster,"_",0))
}

#######################################
# Step 5: Print summary
#######################################

cellsius_print_summary = function(expr_dt_melt){
  cat('--------------------------------------------------------\n',
      'Summary of rare cell types\n',
      '--------------------------------------------------------\n\n')
  for(clust in unique(expr_dt_melt$main_cluster)){
    
    if(!any(expr_dt_melt[main_cluster==clust]$gene_cluster!=0)){
      next
    }
    cat('Main cluster: ', clust,  '\n', '---------------\n')
    subclusts = unique(expr_dt_melt[main_cluster==clust & gene_cluster!=0][order(gene_cluster)]$gene_cluster)
    for(subclust in subclusts){
      
      cat('Subcluster: ', subclust, '\n',
          'Number of cells: ', 
          length(unique(expr_dt_melt[main_cluster==clust & 
                                       sub_cluster == paste(clust,subclust,1,sep="_")]$cell_idx)),
          '\n Marker genes: \n')
      
      print(unique(expr_dt_melt[main_cluster==clust & gene_cluster == subclust][,c("gene_id","symbol","description")]))
      cat('\n\n')
      
    }
  }
}

############################################
# OPTIONAL: Final assignment to unique clusters
# Note: This is different from the previous subcluster asignemnt, where a cell can potentially be
# a member of multiple subgroups.
############################################

cellsius_final_cluster_assignment = function(rare, sce, group_id, min_n_genes = 3){
  
  rare[,n_genes:=length(unique(gene_id)),by='sub_cluster']
  
  assignments = data.table(cell_idx = colnames(sce), pData(sce)[,group_id])
  names(assignments) = c('cell_idx', 'group')
  assignments$group = as.character(assignments$group)
  assignments = merge(assignments, rare[n_genes>=min_n_genes,c('cell_idx','main_cluster','sub_cluster')],by='cell_idx',all=T)
  assignments = unique(assignments)
  
  final_assignment = function(main_cluster,sub_cluster){
    
    if(length(sub_cluster)==1){
      if(is.na(sub_cluster) || grepl("0$",sub_cluster)){
        out = main_cluster
      } else {
        out = gsub('_\\d$','',sub_cluster)
      }
    } else {
      subclusts = gsub('_\\d$', '',sub_cluster[grepl("1$",sub_cluster)])
      out = paste(subclusts,collapse='-')
      if(out == ''){out = main_cluster}
    }
    return(out)
  }
  
  assignments[,final:=final_assignment(group,sub_cluster),by="cell_idx"]
  assignments = unique(assignments[,c('cell_idx','final')])
  
  out = data.frame(cluster = as.character(assignments$final), row.names = assignments$cell_idx)
  out = out[colnames(sce),,drop=F]

  return(out)
}


# Visualize output of CellSIUS on tSNE map
# tsne = RTsne object 
# rare = output of the rare cell types algorithm (a data.table)
plot_rare_cells = function(tsne,rare){
  
  tsne_dt = data.table(tSNE1 = tsne$Y[,1], tSNE2 = tsne$Y[,2], cell_idx = rownames(tsne$Y))
  tsne_dt = merge(tsne_dt, rare[,c('cell_idx','main_cluster','sub_cluster')],
                  by = c('cell_idx'), all = T)
  tsne_dt[is.na(main_cluster),main_cluster:='Other']
  tsne_dt[main_cluster == 'Other',sub_cluster:='none']
  tsne_dt[grepl('_0$',sub_cluster),sub_cluster:= 'none']
  
  setkey(tsne_dt, 'cell_idx')
  tsne_dt = unique(tsne_dt)
  
  rc_cols = brewer.pal(10,"Spectral")[rep(c(1,9,7,2,6,10,3,8),3)]
  
  p = ggplot(tsne_dt, aes(x = tSNE1, y= tSNE2)) + 
    geom_point(color = "darkgray", alpha = 0.5, size = 1.5)+
    theme_bw() + theme(text = element_text(size = 15))
  p = p + geom_point(data = tsne_dt[sub_cluster!='none'], aes(x=tSNE1, y=tSNE2, color = sub_cluster))+
    scale_color_manual(values = rc_cols) + guides(color = guide_legend(title = 'Subcluster'))
  return(p)                                                                     
}

# legacy support
rare_cell_type_identifier = function(sce,group_id,min_n_cells=10,verbose = T, min_fc = 2,
                                     organism = "human", corr_cutoff = NULL, iter=0, max_perc_cells = 50,
                                     fc_between_cutoff = 1){
  return(cellsius_main(sce, group_id, min_n_cells, verbose, min_fc, organism, corr_cutoff, iter,
                       max_perc_cells, fc_between_cutoff))
}
