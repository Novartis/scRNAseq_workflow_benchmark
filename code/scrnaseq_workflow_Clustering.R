#####################################################
#
# scRNASeq pipeline functions
#
# PART V: Clustering
# _______________________
#
# This script contains wrapper functions to different clustering methods.
# 
# Authors:
#   Rebekka Wegmann (rebekka.wegmann@novartis.com)
#   Marilisa Neri (marilisa.neri@novartis.com) 
####################################################




#######################################
# Clustering
#######################################

dist.gen =  function(x,method="euclidean", ...) if ( method %in% c("spearman","pearson","kendall") ) as.dist( (1 - cor(t(x),method=method,...))/2 ) else dist(x,method=method,...)

#________________________________
# MCL 
# Note that this needs the mcl binary installed externally
# MCL can be found here:http://micans.org/mcl/
# Input:
# - mat = matrix of normalized expression values on log scale (i.e. norm_exprs slot
#   of the SCESet) or a similarity matrix, id is_similarity = TRUE
#_______________________________

build_adjacency_matrix = function(mat,cutoff="auto", is_similarity = F){
  library(Ckmeans.1d.dp)
  if(!is_similarity){
    message("Computing cell Pearson correlation coefficient")
    corr.cells = cor(mat,method="pearson")
  } else {corr.cells = mat}
  
  adj.corr = corr.cells
  
  if(cutoff=="auto"){
    # we find the best correlation cutoff by looking for a "valley"
    # in the histogram of correlations. This function attempts to set the
    # cutoff automatically, but might not always succeed...
    
    # if there are more than 500 cells, randomly sample 500 correlations
    if(dim(corr.cells)[1]>500){
      idx = sample(seq(dim(corr.cells)[1]),size=500)
    } else {idx = seq(dim(corr.cells)[1])}
    
    freqs = hist(corr.cells[idx,idx],breaks=dim(corr.cells[idx,idx])[1]/10)
    k1d = Ckmeans.1d.dp(corr.cells,k=2)
    cutoff = max(as.vector(corr.cells)[which(k1d$cluster==1)])
    abline(v=cutoff,col="red")
  } else if (is.numeric(cutoff)){cutoff=cutoff} else {
    stop("Please provide a numeric value for corr.cutoff or set to \"auto\"")
  }
  
  message("Building the adjacency matrix")
  adj.corr[adj.corr<cutoff]=0
  adj.corr[adj.corr>0] = 1
  return(list(adj=adj.corr,cor=corr.cells,cutoff=cutoff))
}

MCLcell.clust=function(adj_list,selfloop=T,mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl"){
  library(igraph)
  
  adj = adj_list$adj
  corr.cells = adj_list$cor
  corr.cutoff = adj_list$cutoff
  
  if(!selfloop) diag(adj)=0 # no self-loop for MCL
  message("Building Graph")
  graphs = get.data.frame( graph.adjacency(adj), what = "edges") # gene connection for graphs
  graphs = data.frame(graphs,CORR=sapply(seq(dim(graphs)[1]), function(i) corr.cells[graphs$from[i],graphs$to[i]] -corr.cutoff))
  
  write.table(graphs, file = "tmp.mcl.inp",row.names=F,col.names=F,sep = " ")
  message("Running MCL")
  system(paste0(mcl_path, " tmp.mcl.inp --abc -o tmp.mcl.out"))
  x = scan("tmp.mcl.out", what="", sep="\n")
  MCL.cells = strsplit(x, "[[:space:]]+")
  MCL.cells = lapply(seq(length(MCL.cells)), function(i){
    tmp = sapply(seq(length(MCL.cells[[i]])),function(j){
      gsub('\"','',MCL.cells[[i]][j])
    })
  })
  system("rm tmp.mcl.inp tmp.mcl.out")
  
  groups.MCL = matrix(rep(-1,dim(corr.cells)[2]),ncol=1)
  rownames(groups.MCL) = colnames(corr.cells)
  for(i in seq(length(MCL.cells))) groups.MCL[MCL.cells[[i]],]=i
  
  #if necessary, collapse all clusters containing only 1 cell to a big "unassigned"
  groups.MCL[groups.MCL %in% names(table(groups.MCL)[which(table(groups.MCL)==1)])] = 0
  
  return(groups.MCL)
}

#___________________________________________________
# DBSCAN
# Note that this has issues in high-dimensional space
# You should therefore use as input"
# - dist = the cell-to-cell distance in PCA space
# - min_pts should be set to number of used PCs +1 (don't use more than ~10 PCs)
# - eps is the y value where the elbow on the knn_plot is. Set it to "auto" if you
#   want the function to automatically determine eps
#_________________________________________________
run_dbscan = function(dist,eps="auto",min_pts,tol=0.01){
  library(dbscan)
  #automatic determination of eps (the "elbow" in the kNNdistplot)
  if(eps=="auto"){
    kNNdist = sort(kNNdist(dist,min_pts))
    i = seq(1,length(kNNdist),as.integer(0.001*length(kNNdist)))
    slope_prev = 100
    for(indx in i){
      slope = kNNdist[indx]/indx
      if(slope_prev>=slope-tol*slope){
        slope_prev = slope
      } else {
        elbow = indx
        break
      }
    }
    eps = kNNdist[elbow]
    print(paste("Epsilon: ",eps))
  } else if(!is.numeric(eps)){
    stop("Please provide a value for eps or set it to \"auto\"")} else {eps=eps}
  
  kNNdistplot(dist,k=min_pts)
  abline(h=eps,col="red")
  res = dbscan(dist,eps = eps,minPts = min_pts)
  return(res$cluster)
}

#__________________________________________
# Mclust
#__________________________________________

gmm_main = function(norm_exprs=NULL,pc_scores=NULL,n_comp=10,do_crossval=T, model_type = "VVI",
                     best_k=NULL,tolerance = 1.96,k=1:10,n_cores = 4, return_model = F){
  library(mclust)
  library(parallel)
  
  if(is.null(pc_scores)){
    if(is.null(norm_exprs)){
      stop("Missing expression values. Please provide either a matrix of normalized counts or pre-computed PCA scores.")
    }
    print("Running PCA...")
    pca = pca_stuff(norm_exprs)
    top_pc_scores = pca$x[,1:n_comp]
  } else {
    top_pc_scores = pc_scores[,1:n_comp]
  }
  
  if(do_crossval){
    #fit models with k-fold crossvalidation
    folds = 1:10
    n_folds = length(folds)
    # this randomly determines which samples should be excluded from
    # model fitting during each fold of the cross-validation. 
    idx = 	sample(rep(1:n_folds, length.out = nrow(top_pc_scores))) 
    
    #Set up parallel processing
    library(parallel)
    cl = makeCluster(n_cores) 
    funs = as.character(lsf.str(envir=parent.frame())) #let clusters know about functions in workspace
    clusterExport(cl,funs)
    
    #benchmark
    #time = system.time(parSapply(cl,folds,cross_val,data=testset,idx=idx,structure='VVI',components=k))
    print("Determining number of clusters...")
    likes = parSapply(cl,folds,cross_val,data=top_pc_scores,idx=idx,structure=model_type,components=k)
    stopCluster(cl)
    
    mean_likes = apply(likes,1,function(x) sum(x[which(is.finite(x))])/length(which(x!=0 & is.finite(x))))
    sd_likes = apply(likes,1,function(x) sd(x[which(x!=0 & is.finite(x))]))
    sd_likes = sd_likes[which(!is.na(sd_likes))]
    mean_likes = mean_likes[which(!is.na(sd_likes))]
    best_idx = which(mean_likes==max(mean_likes))
    ci_likes = mean_likes[best_idx]-tolerance*sd_likes[best_idx]
    
    best_k_idx = min(which(mean_likes>=ci_likes)) #smallest numebr of components that
    #fit reasonably well
    best_k = k[best_k_idx]
    mean_likes[best_k_idx]
    col=rep(1,n_folds)
    col[best_k_idx] = 33
    
    # Plot likelihood vs number of clusters
    # This should look like a ROC curve ideally. The best model should have
    # a high likelihood with the smallest no. of clusters, i.e. be the one
    # where the slope decreases. If there is no clear decrease in the
    # slope, this means that pretty much any model fits equally well/bad and that
    # most likely the clustering produces a nonsensical result
    like_dt = data.table(cluster=k,color=as.factor(col),likes)
    like_dt_melt = melt(like_dt,id.vars=c("cluster","color"),val="log-Likelihood",var="fold")
    p = ggplot(like_dt_melt[`log-Likelihood`!=0 & is.finite(`log-Likelihood`)],aes(x=as.factor(cluster),y=`log-Likelihood`,fill=color)) +
      geom_boxplot() + labs(x = "Number of clusters",
                            title = paste0("No. clusters v.s. log-likelihood, ",n_folds,"-fold crossvalidation"),
                            subtitle = paste0(n_comp," principal components used to calcualte model")) + 
      theme_bw() +scale_fill_manual(values=c("white","red"),guide=F)
    
    #ggsave(p,file=file.path(plotdir,paste0("mclust_crossval_",n_comp,".pdf")),height=5,width=7)
    print(p)
    
    print(paste0("Found ",best_k," clusters."))
  } else {best_k = best_k}
  
  if(is.null(best_k)){
    stop("Please provide a value for best_k or set do_crossval = TRUE")
  }
  print("Assigning cells to clusters...")
  #assignments of cells to clusters
  model = calc_model(top_pc_scores,best_k,model_type)
  cluster_assignments = model$classification
  if(return_model){
    return(model)
  } else {
    return(cluster_assignments)}
}

##################################################################################
#Functions called by gmm_main
##################################################################################
#setup paths and stuff
#do pca
pca_stuff = function(log_data_hv,scale_pca=T,center_pca=T){
  pca = prcomp(t(log_data_hv[,-1,with=F]),scale=scale_pca,center=center_pca)
  return(pca)
}

#Fitting GMM with mclust:
##############################################################################
# function to calculate different models
# k = number of compnents
# structure = model structure / constraints. See mclustModelNames for details.

calc_model = function(data,k,structure){
  return(Mclust(data,G=k,modelNames = structure,initialization=
                  list(subset=sample(1:nrow(data),size=as.integer(nrow(data)*4/5)))))
}

#############################################################################
#Functions to calculate log-likelihood out of what mclust returns

# Probability density function for a Gaussian mixture
# Presumes the mixture object has the structure used by mclust

dnormalmix = function(x,mixture,log=FALSE) {
  lambda 	= mixture$parameters$pro
  k 		= length(lambda)
  # Calculate share of likelihood for all data for one component
  like_component = function(x, component) {
    lambda[component] * dmvnorm(
      x,
      mean = mixture$parameters$mean[,component],
      sigma 	= mixture$parameters$variance$sigma[,,component]
    )
  }
  # Create array with likelihood shares from all components over all data
  likes 	= sapply(1:k, like_component ,x = x)
  # Add up contributions from components
  d 	= rowSums(likes)
  if (log) {
    d 	= log(d)
  }
  return(d)
}

# Log likelihood function for a Gaussian mixture, potentially on new data
loglike_normalmix = function(x,mixture) {
  loglike 	= dnormalmix(x, mixture, log = TRUE)
  return(sum(loglike))
}

###############################################################################
#Cross validation things
#Cross validation
#data = input data
#idx = a random sample of folds (e.g. 11432...)
#fold = the current fold
#structure = model structure for mclust (e.g. 'VVI)
#components = a vector containing the numebr of components for which to test models

cross_val = function(fold,data,idx,structure,components){
  #library(mclust,lib.loc = .libPaths()[[2]]) #for the VM
  library(mclust)
  library(mvtnorm)
  like_test = c()
  for(k in components){
    out = tryCatch(
      {
        calc_model(data[which(idx!=fold),],k,structure)
      },
      error=function(cond){
        #try to find another model
        out2 = 0
        counter = 0
        while(out2 == 0 && counter<=5){
          out2 = tryCatch(
            {
              calc_model(data[which(idx!=fold),],k,structure)
            },
            error = return(0),
            warning=return(0),
            finally= {counter = counter +1}
          )
        }
        message('There was an error: \n')
        message(cond)
        write.csv(cond,'gmm.log',append=TRUE)
        return(out2)
      },
      warning=function(cond){
        #try to find another model
        out2 = 0
        counter = 0
        while(out2 == 0 && counter<=10){
          out2 = tryCatch(
            {
              calc_model(data[which(idx!=fold),],k,structure)
            },
            error = return(0),
            warning=return(0),
            finally={counter = counter +1}
          )
        }
        message('There was a warning: \n')
        message(cond)
        write.csv(cond,'gmm.log',append=TRUE)
        return(out2)
      },
      finally={
        message('\n done.')
      }
    )
    if(class(out)=='Mclust'){  
      like_test = append(like_test,loglike_normalmix(data[which(idx==fold),],out))
    }
    else{
      like_test = append(like_test,0)
    }
  }
  return(like_test)
}

#___________________________________________
# SEURAT
#___________________________________________
seurat_clustering = function(sce,vars.to.regress=NULL,res=0.6,n_comp=10){
  
  library(Seurat)
  #make SEURAT object, scale and optionally regress out confounders
  tmp_seurat = CreateSeuratObject(raw.data = counts(sce))
  tmp_seurat@data = norm_exprs(sce) #add the normalized values
  
  # This next step is a bit of cheating. Seurat expects us to run the complete
  # workflow on the same object and checks whether data have been normalized
  # by checking if object@calc.params$NormalizeData$normalization.method exists.
  # Since we provided normalized values, and do not want to re-run normalization,
  # we just put a dummy value in that slot.
  
  tmp_seurat@calc.params$NormalizeData = list(normalization.method ="dummy")
  
  if(!is.null(vars.to.regress)){
    if(any(!vars.to.regress%in%names(pData(sce)))){
      stop("Variables to regress out have to be column names in pData(sce)")
    }
    tmp_seurat = AddMetadata(object = tmp_seurat, metadata = pData(sce)[,vars.to.regress])
    tmp_seurat = ScaleData(object = tmp_seurat,vars.to.regress=vars.to.regress)
  } else {
    tmp_seurat = ScaleData(object = tmp_seurat)
  }
  
  tmp_seurat = RunPCA(object = tmp_seurat, pc.genes = rownames(sce), do.print = FALSE) 
  tmp_seurat = FindClusters(object = tmp_seurat, reduction.type = "pca", dims.use = 1:n_comp, 
                            resolution = res, print.output = 0, save.SNN = TRUE)
  seurat_assignment = tmp_seurat@ident
  return(seurat_assignment)
}


