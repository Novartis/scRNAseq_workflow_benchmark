##################################################################################################
# scRNAseq workflow - main file
#
# This file contains all the code (without explanation) that is displayed in the html vignette.
# It might be easier for you to run the workflow from the .Rmd file in the vignettes folder,
# because this contains both the code and extensive documentation.
#
# Author: Rebekka Wegmann and Marilisa neri
#############################################################################################

#-----------------------------------------------
#                 LICENSE
#-----------------------------------------------
# Copyright 2018 Novartis Institutes for BioMedical Research Inc.
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



## ------------------------------------------------------------------------
#source("https://bioconductor.org/biocLite.R")

## ---- eval = F-----------------------------------------------------------
## install.packages("devtools")

## ---- eval=F-------------------------------------------------------------
## # from CRAN
## install.packages("Rtsne") #Rtsne v. 0.13
## install.packages("ggplot2") # ggplot2 v. 2.2.1
## install.packages("data.table") # data.table v. 1.10.4
## install.packages("RColorBrewer") # RColorBrewer v. 1.1-2
## install.packages("mvoutlier") # mvoutlier 2.0.8, required by some functions from scater
## 
## devtools::install_github("bwlewis/irlba") # irlba 2.3.2, optional. Make sure you use irlba > 2.3.1, older versions contain a bug that results in unreliable output!
## 
## # rom Bioconductor
## biocLite("scater") # scater v. 1.4.0
## biocLite("scran") # scran v. 1.4.5

## ---- eval=F-------------------------------------------------------------
## # ensembldb 2.0.4 and EnsDb.Hsapiens.v79 2.1.0
## biocLite(c("ensembldb","EnsDb.Hsapiens.v79"))
## # org.Hs.eg.db v. 3.4.1
## biocLite("org.Hs.eg.db")

## ---- eval=F-------------------------------------------------------------
## # ensembldb 2.0.4 and EnsDb.Mmusculus 2.1.0
## biocLite(c("ensembldb","EnsDb.Mmusculus.v75"))
## # org.Mm.eg.db v. 3.4.1
## biocLite("org.Mm.eg.db")

## ----eval = F------------------------------------------------------------
## # M3Drop version 3.05.00 (Note: This is still under active development, please let me know if a new version breaks my functions...)
## devtools::install_github("tallulandrews/M3D")

## ---- eval = F-----------------------------------------------------------
## install.packages("cluster") # cluster v. 2.0.6
## install.packages("dendextend") # dendextend v. 1.5.2
## install.packages("Ckmeans.1d.dp") # Ckmeans.1d.dp v. 4.2.1
## install.packages("dbscan") # DBSCAN v. 1.1-1
## install.packages(c("mclust","mvtnorm")) # mclust v. 5.4 and mvtnorm 1.0.6
## install.packages("dynamicTreeCut") # dynamicTreeCut v. 1.63-1
## 
## biocLite("SC3") # SC3 v. 1.4.2
## devtools::install_github("satijalab/seurat") #Seurat 2.0.1
## devtools::install_github('JustinaZ/pcaReduce') #pcaReduce 1.68.0

## ----eval = F------------------------------------------------------------
## biocLite("limma") # limma v. 3.32.5
## biocLite("MAST") # MAST v. 1.2.1

## ---- message = FALSE----------------------------------------------------
rm(list=ls())
graphics.off()

wd = ".."

#Directory where input files are stored
input_dir = file.path(wd,"example_data")

#Directory where code is stored
code_dir = file.path(wd,"code")

#where to save the output data?
out_data_dir = file.path(wd,"example_output")
if(!dir.exists(out_data_dir)) {dir.create(out_data_dir)}

#where to save the produced plots?
plotdir = file.path(wd,"example_plots")
if(!dir.exists(plotdir)) {dir.create(plotdir)}

set.seed(17) #to make tSNE plots reproducible

source(file.path(code_dir,"scRNASeq_pipeline_functions.R"))

# loading the libraries that are required throughout the analysis
library_calls()

## ----message=F,eval=F----------------------------------------------------
## # Read the dataset
## counts = read.delim(file.path(input_dir,"pbmc_example_counts.txt"),sep="\t")
## 
## # if we have genes that are not expressed in any cell, discard them
## keep_feature = !gene_filter_by_feature_count(counts,0)
## counts = counts[keep_feature,]
## 
## # make a table of metadata (e.g. batch, cell type annotation,treatment,...)
## # Here, we do not have any such information, so we just give each cell a name
## # Note that the rownames of this table have to correspond to the column names
## # of the count matrix.
## annot = data.frame(cell_idx=paste0("pbmc_C",seq(dim(counts)[2])))
## rownames(annot) = colnames(counts)
## pd = new("AnnotatedDataFrame", data=annot)
## 
## # get gene annotations from ensembldb
## # optional: get gene descriptions from org.Hs.eg.db (slows down the process a lot!)
## # NOTE: the output of get_gene_annotations is a data.table sorted by gene identifier.
## #       This means the genes are no longer in the same order as in the count matrix!
## geneName = get_gene_annotations(rownames(counts),organism = "human",get_descriptions = F)
## 
## # convert this to feature metadata
## fd_table = as.data.frame(geneName)
## rownames(fd_table) = geneName$gene_id
## fd_table = fd_table[rownames(counts),]
## fd = new("AnnotatedDataFrame",data=fd_table)
## 
## #construct SCESet
## sce = newSCESet(countData = counts, phenoData=pd,featureData = fd)
## 

## ----eval=F--------------------------------------------------------------
## #calculate QC metrics
## sce = calculateQCMetrics(sce,feature_controls = list(MT=which(fData(sce)$chr=="MT")))
## 
## # assign cell cycle phase (based on the method from scran)
## # because PBMCs are post-mitotic, most cells should be assigned to G0/G1 phase
## cc = annotate_cell_cycle(sce)
## sce$cell_cycle_phase = cc$phases
## sce$cc_scores = cc$scores
## 
## #save the SCESet
## save(sce,file=file.path(out_data_dir,"sce_raw.RData"))

## ------------------------------------------------------------------------
load(file.path(out_data_dir,"sce_raw.RData"))

p1 = plotQC(sce, type="high",feature_names_to_plot = "symbol") 
print(p1)


## ------------------------------------------------------------------------
p1.2 = custom_plotHighestExprs(sce,feature_names_to_plot = "symbol")
p1.2 = p1.2 + xlab("Expression [raw counts, log2]")
print(p1.2)

# to save a plot, use the ggsave function:
ggsave(p1.2, file = "saved_example_plot.pdf",height=7,width=7)

## ----warning=F-----------------------------------------------------------
p2 = plotQC(sce, type = "exprs-freq-vs-mean") 
p2 = p2+xlab("Mean expression [raw counts, log2]")+ggtitle("Mean expression versus detection rate")
print(p2)
# Check total number of zeroes
t = table(counts(sce)==0)
print(t/sum(t)*100)

## ------------------------------------------------------------------------
p3 = plotQC(sce, type="find", variable="total_features") 
print(p3 + ggtitle("Correlation of principal components with total detected features."))

## ------------------------------------------------------------------------
p3.2 = plotQC(sce, type="find", variable="cell_cycle_phase")
print(p3.2+ggtitle("Correlation of principal components with cell cycle phase."))

## ------------------------------------------------------------------------
vars = c("total_counts","total_features","cell_cycle_phase")
p4 = plotQC(sce, type="expl", variables=vars)
print(p4 + ggtitle("Percentage of explained variance"))

## ------------------------------------------------------------------------
min_genes = 9.0 #minimum number of features (genes) per cell [log2]
min_UMI = 10.5  #minimum total UMIs / cell [log2]
mt_threshold = 9 #Maximum percentage of mitochondrial genes

## ------------------------------------------------------------------------
plot_RNA_QC(sce, min_genes = min_genes, min_UMI = min_UMI)
plot_MT_QC(sce,mt_threshold)


## ------------------------------------------------------------------------
sce$keep_manual = (   !cell_filter_by_feature_count(counts(sce),2^min_genes) &
                      !cell_filter_by_total_UMI(counts(sce),2^min_UMI) &                                !cell_filter_by_mt_content(sce$pct_counts_feature_controls_MT,mt_threshold))

table(sce$keep_manual)

sce_clean = sce[,sce$keep_manual]

## ------------------------------------------------------------------------
n_th = 1
min_counts = 2

keep_feature = !(gene_filter_by_feature_count(counts(sce_clean),n_th, min_counts))
sce_clean = sce_clean[keep_feature,]

## ------------------------------------------------------------------------

sce_clean = calculateQCMetrics(sce_clean,feature_controls = list(MT=which(fData(sce_clean)$chr=="MT")))

# The variables used to detect outliers
vars = c( "pct_counts_top_100_features", 
          "total_features", "pct_counts_feature_controls", 
           "log10_counts_endogenous_features", 
           "log10_counts_feature_controls")

sce_clean = plotPCA(sce_clean,
                    size_by = "total_features", 
                    pca_data_input = "pdata",
                    selected_variables = vars,
                    detect_outliers = TRUE,
                    return_SCESet = TRUE)

table(sce_clean$outlier)

#here, we remove the outliers
sce_clean = sce_clean[,!sce_clean$outlier]

# Finally, we again remove genes that are not expressed
keep_feature = !(gene_filter_by_feature_count(counts(sce_clean),n_th,min_counts))
sce_clean = sce_clean[keep_feature,]

save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))

## ------------------------------------------------------------------------
p = my_plot_PCA(counts = log2(counts(sce_clean)+1),
                scale_pca = T, center_pca = T, return_pca = F, use_irlba=F,
                color = pData(sce_clean)[,"total_features",drop=F])
p = p+ggtitle("PCA on raw log2(counts)")
print(p)

## ------------------------------------------------------------------------
p = my_plot_tSNE(counts = log2(counts(sce_clean)+1),
                 is_distance = F, scale_pca = F, n_comp = 50, return_tsne=F,
                 color = pData(sce_clean)[,"total_features",drop=F])
p = p+ggtitle("t-SNE on raw log2(counts)")
print(p)

## ----message=F-----------------------------------------------------------
# normalize data
sce_clean = normalize_counts(sce_clean,method = "scran")

# normalized values are automatically log2-transformed and
# stored in the norm_exprs slot of the SCESet
norm_exprs(sce_clean)[1:5,1:5]

save(sce_clean, file = file.path(out_data_dir,"sce_clean.RData"))

## ------------------------------------------------------------------------
# PCA of normalized values
sum_norm = data.frame(sum_expression = colSums(norm_exprs(sce_clean)))
p1 = my_plot_PCA(counts = norm_exprs(sce_clean),
                 color=sum_norm,
                 return_pca = F, scale_pca = T, center_pca = T)
p1 = p1+ggtitle("PCA on normalized counts")
print(p1)

#tSNE of normalized values
p2 = my_plot_tSNE(counts = norm_exprs(sce_clean),
                  color=sum_norm,
                  return_tsne = F, is_distance = F)
p2 = p2 + ggtitle("t-SNE (50 PCs) on normalized counts")
print(p2)

## ----message=F-----------------------------------------------------------
cd20 = t(norm_exprs(sce_clean["ENSG00000156738",]))
colnames(cd20) = "CD20"

go_id = "GO:0002376" 
ens_go = GO_to_gene(go_id)
info_GO = rownames(sce_clean)%in%ens_go
table(info_GO)

p = my_plot_PCA(counts = norm_exprs(sce_clean[info_GO,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - GO:0002376 features",
                color = cd20)
print(p)

## ------------------------------------------------------------------------
info_HVG = info.genes(2^norm_exprs(sce_clean)-1,PLOT=T,qcv=0.25,pv=.1,q=.5,minBiolDisp = 0) 
table(info_HVG)
p = my_plot_PCA(counts = norm_exprs(sce_clean[info_HVG,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - HVG features",
                color = cd20)
print(p)

## ------------------------------------------------------------------------
info_NBdrop = run_DANB(counts(sce_clean),method = "NBDrop",save_plot=F, cutoff = 0.1) 
info_NBdisp = run_DANB(counts(sce_clean),method = "NBDisp",save_plot=F, perc_genes = 10) 
table(info_NBdrop,info_NBdisp)

p = my_plot_PCA(counts = norm_exprs(sce_clean[info_NBdrop,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - NBDrop features",
                color = cd20)
print(p)
p = my_plot_PCA(counts = norm_exprs(sce_clean[info_NBdisp,]),
                 return_pca = F, scale_pca = T, center_pca = T,
                title = "PCA - NBDisp features",
                color = cd20)
print(p)

## ---- eval=F-------------------------------------------------------------
## sce_info = sce_clean[info_NBdrop,]
## dim(sce_info)
## 
## # tSNE map of the cleaned data
## # note that by setting return_tsne = T, we can obtain the t-SNE object for later use
## tsne_info = my_plot_tSNE(counts = norm_exprs(sce_info),
##                          scale_pca = F, n_comp = 50, return_tsne=T)$tsne

## ----eval=F--------------------------------------------------------------
## save(sce_info, file = file.path(out_data_dir,"sce_info.RData"))
## save(tsne_info,file = file.path(out_data_dir,"tsne_info.RData"))

## ------------------------------------------------------------------------
#load the data we need
load(file.path(out_data_dir,"sce_info.RData"))
load(file.path(out_data_dir,"tsne_info.RData"))

## ---- fig.align='center', fig.width=12, fig.height=8---------------------
b_cell = t(norm_exprs(sce_info["ENSG00000156738",]))
colnames(b_cell) = "B-cell"

monocyte = data.frame(Monocyte = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('CD14','LYZ','FCGR3A','MS4A7')),]))

t_cell = data.frame(`T-cell` = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('CD3E','CD3D','CD3G')),]))

nk_cell = data.frame(`NK cell` = colSums(norm_exprs(sce_info)[which(fData(sce_info)$symbol %in% c('GNLY','NKG7')),]))

# Make plots
# Note that by providing the tsne input variable instead of counts,
# we can use an existing t-SNE calculation for plotting

p1 = my_plot_tSNE(tsne = tsne_info, color = b_cell, title = "B-cell marker expression")  
p2 = my_plot_tSNE(tsne = tsne_info, color = monocyte, title = "Monocyte marker expression")
p3 = my_plot_tSNE(tsne = tsne_info, color = t_cell, title = "T-cell marker expression")
p4 = my_plot_tSNE(tsne = tsne_info, color = nk_cell, title = " NK cell marker expression")
ggmultiplot(p1,p2,p3,p4,cols=2)

## ------------------------------------------------------------------------
assignment = data.table(tsne1 = tsne_info$Y[,1], tsne2 = tsne_info$Y[,2],cell_type = 'T-cell')
assignment[tsne1 < -10 ,cell_type:='B-cell']
assignment[tsne1 > 5 ,cell_type:='Monocyte']
assignment[tsne2 < -17 & tsne1 > -1,cell_type:='NK Cell']

sce_info$cell_type = assignment$cell_type

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"cell_type",drop=F])
print(p+labs(title="t-SNE on informative genes",subtitle = "Colored by manual cell annotation"))


## ---- echo=F-------------------------------------------------------------
library(SC3)

## ----eval=F--------------------------------------------------------------
## library(SC3)
## sce_info = sc3_prepare(sce_info, ks = 2:10, n_cores = 4)

## ----eval=F--------------------------------------------------------------
## sce_info = sc3_estimate_k(sce_info)
## sce_info@sc3$k_estimation

## ----eval=F--------------------------------------------------------------
## sce_info = sc3(sce_info, ks = 6, biology = TRUE, n_cores = 4)

## ---- eval = F-----------------------------------------------------------
## sce_info = sc3(sce_info, ks = 6, biology = F, n_cores = 8)
## sce_info = sc3_run_svm(sce_info)
## 
## sce_info@sc3$svm_train_inds = NULL
## sce_info = sc3_calc_biology(sce_info, k=c(8,13), regime = "marker")
## 
## # to visualize the markers, use my modified fuction:
## 
## # change the plotted gene names to symbol for better readability
## plot_sce = sce_info
## rownames(plot_sce) = fData(plot_sce)$symbol
## custom_sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90)
## 
## rm(plot_sce)

## ---- out.width="110%",fig.height=4--------------------------------------
sc3_plot_consensus(sce_info, k=6)

## ---- out.width="110%",fig.height=4--------------------------------------
sc3_plot_expression(sce_info, k = 6, show_pdata = c("cell_type"))

## ---- out.width="110%",fig.height=6--------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info
rownames(plot_sce) = fData(plot_sce)$symbol
sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90, show_pdata = c("cell_type"))

# for hybrid SVM approach:
# custom_sc3_plot_markers(plot_sce, k=6, p.val = 0.01, auroc = 0.90)

## ------------------------------------------------------------------------
assignment = data.table(clust = sce_info$sc3_6_clusters, cell_type = sce_info$cell_type)
assignment[clust == 1, cell_type:= 'T helper cell']
assignment[clust == 2, cell_type:= 'Monocyte']
assignment[clust==3,cell_type:='?']
assignment[clust==4,cell_type:='CD8+ T-cell']
assignment[clust==5,cell_type:='NK cell']
assignment[clust == 6, cell_type:= 'B-cell']

sce_info$SC3_assignment = assignment$cell_type
p = my_plot_tSNE(tsne = tsne_info,
                 color = pData(sce_info)[,"SC3_assignment",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "SC3 assignment",
                 show_proportions = T)
print(p)
save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))

## ------------------------------------------------------------------------
pca = my_plot_PCA(counts = norm_exprs(sce_info),return_pca=T)$pca
screeplot(pca,type = 'lines')

## ------------------------------------------------------------------------
dist_eucl = dist.gen(pca$x[,1:6],method='euclidean')
hfit = hclust(dist_eucl,method="average")
plot(hfit, labels = F, main = "hclust on euclidean distance in PCA space")

## ---- message =F---------------------------------------------------------
library(dynamicTreeCut)
groups_hclust_eucl = cutreeDynamic(hfit, distM = as.matrix(dist_eucl), deepSplit = 0, minClusterSize = 5, maxCoreScatter = 0.70, minGap = 0.25)

## ------------------------------------------------------------------------
library(cluster)
si = silhouette(groups_hclust_eucl,dist_eucl)
plot(si, col = "darkgray", border=NA, main = "Silhouette plot for hclust (euclidean in PCA space)")

## ------------------------------------------------------------------------
sce_info$hclust_sil = si[,3]
p = my_plot_tSNE(tsne = tsne_info,
                 color = pData(sce_info)[,"hclust_sil",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "tSNE colored by silhouette width")
print(p+scale_color_distiller(type="div", palette = "RdBu"))

## ------------------------------------------------------------------------
sce_info$hclust_eucl = as.factor(groups_hclust_eucl)
table(sce_info$SC3_assignment,sce_info$hclust_eucl)

## ------------------------------------------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_info
rownames(plot_sce) = fData(plot_sce)$symbol
monocyte_markers = which(fData(sce_info)$symbol %in% c('CD14','LYZ','FCGR3A','MS4A7'))
p = plotExpression(plot_sce,features = monocyte_markers, x="hclust_eucl",
               colour_by = "hclust_eucl")
print(p)

## ----fig.width = 8-------------------------------------------------------
sce_info$hclust_eucl = factor(sce_info$hclust_eucl, levels = levels(sce_info$hclust_eucl), labels = c("T-helper cell","CD14++/CD16- Monocyte","B-cell","CD8+ T-cell","NK cell","CD14+/CD16++ Monocyte"))
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"hclust_eucl",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (euclidean) assignment",
                 show_proportions = T)
print(p)

## ------------------------------------------------------------------------
library(dendextend)
dist_pearson = dist.gen(t(norm_exprs(sce_info)),method = "pearson")
hfit = hclust(dist_pearson,method="average")
groups_hclust_pearson = cutree(hfit, k=6)
sce_info$hclust_pearson = as.factor(groups_hclust_pearson)
hfit2 = color_branches(hfit, k=6)
hfit2 = hfit2 %>% set("labels", rep("",dim(sce_info)[2]))
plot(hfit2, main = "hclust on Pearson correlation")

## ----fig.width=8---------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"hclust_pearson",drop=F],
                 shape = pData(sce_info)[,"cell_type",drop=F],
                 title = "hclust (Pearson) assignment",
                 show_proportions = T)
print(p)

## ------------------------------------------------------------------------
library(pcaReduce)
expr_mat_info = t(norm_exprs(sce_info))
assignment_info = PCAreduce(expr_mat_info, nbt = 5, q = 6, method = "M")

## ----fig.width=12,fig.height=6-------------------------------------------
table(assignment_info[[1]][,4])

plots = list()
for(i in 1:5){
  df = data.frame(apply(assignment_info[[i]],c(1,2),as.factor))
  p = my_plot_tSNE(tsne=tsne_info, color = df[,'cl_id.2',drop=F],
                   shape = pData(sce_info)[,"cell_type",drop=F],
                   title = paste0("pcaReduce, run ", i))
  plots[[i]] = p
}

ggmultiplot(plots[[1]],plots[[2]],plots[[3]],plots[[4]],plots[[5]],cols=3)


## ---- message = F--------------------------------------------------------
sim = 1-as.matrix(dist_eucl)/max(dist_eucl)
adj_corr = build_adjacency_matrix(mat = sim,cutoff="auto",is_similarity=T)

## ---- message = F--------------------------------------------------------
groups_MCL = MCLcell.clust(adj_corr,mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl")
sce_info$MCL = as.factor(groups_MCL)
table(sce_info$MCL,sce_info$hclust_eucl)

## ---- fig.width = 8------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"MCL",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=auto) assignment",
                 show_proportions = T)
print(p)

## ----fig.height=6,fig.width = 9------------------------------------------
adj_corr = build_adjacency_matrix(mat = sim,cutoff=0.8,is_similarity=T)
groups_MCL = MCLcell.clust(adj_corr,mcl_path = "/da/dmp/cb/prog/mcl-14-137/bin/mcl")
sce_info$MCL2 = as.factor(groups_MCL)
table(sce_info$MCL2,sce_info$hclust_eucl)

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"MCL2",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "MCL (cutoff=0.7) assignment",
                 show_proportions = T)
p = p + scale_color_manual(values = brewer.pal(10,"Paired"))
print(p)

## ----eval = F------------------------------------------------------------
## library(Seurat)
## seurat_assignments = seurat_clustering(sce_info,vars.to.regress=NULL,res=1.2)
## sce_info$seurat = as.factor(seurat_assignments)

## ----fig.width=8, eval = F-----------------------------------------------
## table(sce_info$seurat,sce_info$hclust_eucl)
## p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"seurat",drop=F],
##                  shape = pData(sce_info)[,"hclust_eucl",drop=F],
##                  title = "Seurat assignment",show_proportions = T)
## print(p)

## ------------------------------------------------------------------------
DBSCAN_groups = run_dbscan(dist_eucl,eps="auto",min_pts=7,tol=0.005)

## ----fig.height=6,fig.width=8--------------------------------------------
sce_info$DBSCAN = as.factor(DBSCAN_groups)
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"DBSCAN",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "DBSCAN assignment",show_proportions = T)
print(p)


## ----fig.width=7,fig.height=6--------------------------------------------
si = silhouette(DBSCAN_groups, dist_eucl)
plot(si, col = "darkgray", border=NA, main = "Silhouette plot for dbscan (euclidean in PCA space)")

## ------------------------------------------------------------------------
boot = fpc::clusterboot(data = dist_eucl, clustermethod = fpc::dbscanCBI, eps = 5.53,
                         MinPts = 7,B = 100, distances = T, count=F, method = "dist",
                        bootmethod = "boot")

dt = melt(data.table(cluster = c(0:4), boot$bootresult),id.vars="cluster",val="jaccard index",var="boot number")
dt = dt[cluster!=0]
p = ggplot(dt, aes(x=as.factor(cluster),y=`jaccard index`)) + geom_boxplot() +theme_bw()
print(p + ggtitle("Jaccard similarity between cluster assignment on the full data and on 100 bootstrap samples"))

## ------------------------------------------------------------------------
library(mclust)
mclust_assignment = gmm_main(pc_scores = pca$x,n_comp=6, k=1:8,n_cores=1)

## ------------------------------------------------------------------------
sce_info$mclust = as.factor(mclust_assignment)
table(sce_info$mclust,sce_info$hclust_eucl)

## ---- fig.width = 8------------------------------------------------------
p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"mclust",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "Mclust assignment",show_proportions = T)
print(p)

## ---- fig.width=8--------------------------------------------------------
mclust_model = gmm_main(pc_scores = pca$x,n_comp=6,do_crossval = F, best_k = 6, return_model = TRUE)
sce_info$mclust_forced = as.factor(mclust_model$classification)

table(sce_info$mclust_forced,sce_info$hclust_eucl)

p = my_plot_tSNE(tsne = tsne_info, color = pData(sce_info)[,"mclust_forced",drop=F],
                 shape = pData(sce_info)[,"hclust_eucl",drop=F],
                 title = "Mclust assignment")
print(p)

## ----message=F,results="hide"--------------------------------------------
mclust_boot = MclustBootstrap(mclust_model, nboot=500)

## ------------------------------------------------------------------------
summary(mclust_boot, what = "ci")

## ----fig.width=8, fig.height = 4-----------------------------------------
par(mfrow=c(1,6))
plot(mclust_boot, what = "pro")

## ---- fig.width=8,fig.height=12------------------------------------------
par(mfrow=c(6,6))
plot(mclust_boot, what = "mean")
par(mfrow=c(6,6))
plot(mclust_boot, what = "var")
par(mfrow=c(1,1))

## ---- eval = F-----------------------------------------------------------
## save(sce_info,file = file.path(out_data_dir,"sce_info.RData"))

## ---- warning=F----------------------------------------------------------
sce_clean$cell_type = sce_info$cell_type
cellsius_out = cellsius_main(sce_clean, group_id = "cell_type", min_n_cells = 5,
                                 verbose =T, min_fc = 2, organism = "human", iter=0, 
                                 max_perc_cells = 50, fc_between_cutoff = 1)

## ------------------------------------------------------------------------
cellsius_out

## ------------------------------------------------------------------------
unique(cellsius_out[cell_idx=="AAAGCAAGTCAAGCGA"]$sub_cluster)

## ------------------------------------------------------------------------
unique(cellsius_out[sub_cluster == "T-cell_1_1"][,c('gene_id','symbol','description')])

## ------------------------------------------------------------------------
cellsius_out[,length(unique(cell_idx)),by="sub_cluster"]

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out, tsne=tsne_info)

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out[main_cluster=="T-cell"], tsne=tsne_info)

## ------------------------------------------------------------------------
plot_rare_cells(rare=cellsius_out[sub_cluster=="T-cell_1_1"], tsne=tsne_info)

## ---- fig.width=8,fig.height=6-------------------------------------------
marker_idx = which(rownames(sce_clean)%in%cellsius_out[sub_cluster=='T-cell_1_1']$gene_id)

plotlist = list()
colors = t(norm_exprs(sce_clean)[marker_idx,,drop=F])
colnames(colors) =  fData(sce_clean)$symbol[marker_idx]

for(i in seq(dim(colors)[2])){
  plotlist[[i]] = my_plot_tSNE(tsne = tsne_info,
                      color = colors[,i,drop=F],
                      alpha = 0.8, title = colnames(colors)[i]) + scale_color_distiller(palette = 'RdYlBu')
}

ggmultiplot(plotlist[[1]],plotlist[[2]],plotlist[[3]], plotlist[[4]],cols = 2)

## ---- eval = F-----------------------------------------------------------
## library(shiny)
## launch_marker_vis_app(tsne = tsne_info, sce=sce_clean, marker_idx = marker_idx)

## ------------------------------------------------------------------------
final_clusters = cellsius_final_cluster_assignment(cellsius_out, sce_clean, group_id = "cell_type", min_n_genes=3)

table(final_clusters$cluster)

p = my_plot_tSNE(tsne = tsne_info, color = final_clusters, title = "CellSIUS cluster assignment")
print(p)

## ------------------------------------------------------------------------
sce_clean$SC3_assignment = sce_info$SC3_assignment
sce_clean$hclust_eucl = sce_info$hclust_eucl

## ------------------------------------------------------------------------
de_wilcox = run_wilcoxon_test(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                           fc_cutoff = 1, alpha = 0.05)
dim(de_wilcox)
head(de_wilcox)

## ------------------------------------------------------------------------
table(de_wilcox$DE_flag)

## ------------------------------------------------------------------------
mean_exprs_dt = data.table(gene_id = rownames(sce_clean),mean_exprs = fData(sce_clean)$mean_exprs)
de_wilcox = merge(de_wilcox,mean_exprs_dt, by = "gene_id")
names(de_wilcox)

## ------------------------------------------------------------------------
p = generic_scatterplot(de_wilcox, x_col = "mean_exprs", y_col = "log2fc",
                        color = "DE_flag")
print(p+ggtitle('MA plot for Wilcoxon test'))

## ------------------------------------------------------------------------
de_wilcox[,log10_pval:=-log10(adj_pval)]
p = generic_scatterplot(de_wilcox, x_col = "log2fc", y_col = "log10_pval",
                        color = "DE_flag")
print(p+ggtitle('Volcano plot for Wilcoxon test'))

## ----message=F-----------------------------------------------------------
library(DT)
top_DE = de_wilcox[DE_flag==TRUE][order(log2fc,decreasing = T)]$gene_id[1:20]
gene_table = get_gene_annotations(top_DE,get_descriptions = T)
datatable(gene_table, caption = "Top 20 upregulated genes in B-cells")

## ------------------------------------------------------------------------
de_limma_voom = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "voom", fc_cutoff = 1, alpha = 0.05, count_thr = 1,pct=50)

## ------------------------------------------------------------------------
dim(de_limma_voom)
names(de_limma_voom)

## ------------------------------------------------------------------------
table(de_limma_voom$DE_flag)
table(de_wilcox[de_limma_voom$gene_id,on="gene_id"]$DE_flag,de_limma_voom$DE_flag)

## ------------------------------------------------------------------------
de_limma_voom = merge(de_limma_voom,mean_exprs_dt, by = "gene_id")
names(de_limma_voom)

p = generic_scatterplot(de_limma_voom, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag")
print(p+ggtitle('MA plot for limma-voom'))

## ------------------------------------------------------------------------

voom_no_filt = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "voom", fc_cutoff = 1, alpha = 0.05, count_thr = 0)
voom_no_filt = merge(voom_no_filt, mean_exprs_dt)
p = generic_scatterplot(voom_no_filt, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag")

print(p+ggtitle("MA plot for limma-voom, no filtering"))

table(voom_no_filt[DE_flag==TRUE]$gene_id%in%de_limma_voom[DE_flag==TRUE]$gene_id)

## ------------------------------------------------------------------------
example_dt = copy(voom_no_filt)
example_dt[,logFC:=logFC-1*1/(mean_exprs+1)]
example_dt[,DE_flag := adj.P.Val < 0.05 & abs(logFC)>1]
p = generic_scatterplot(example_dt, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag") + geom_hline(yintercept = 0)
print(p+ggtitle("An example of an MA plot for limma gone wrong"))


## ------------------------------------------------------------------------
de_limma_trend = run_limma(sce_clean, cl_id = "hclust_eucl", cl_ref = "B-cell",
                          method = "trend", fc_cutoff = 1, alpha = 0.05, count_thr = 1, pct=50)

de_limma_trend = merge(de_limma_trend, mean_exprs_dt)
de_limma_trend[,voom_overlap:=de_limma_trend$DE_flag == de_limma_voom$DE_flag]

p = generic_scatterplot(de_limma_trend, x_col = "mean_exprs", y_col = "logFC",
                        color = "DE_flag", shape = "voom_overlap")

print(p + ggtitle("MA plot for limma-trend"))

table(de_limma_trend$DE_flag,de_limma_voom$DE_flag)
table(de_wilcox[de_limma_trend$gene_id,on="gene_id"]$DE_flag,de_limma_trend$DE_flag)

# Compare p-values between Wilcoxon test and limma
ggplot(data.table(wilcox = de_wilcox[de_limma_trend$gene_id,on="gene_id"]$pval,
                  limmatrend = de_limma_trend$P.Val),
       aes(x=-log10(wilcox),y=-log10(limmatrend)))+geom_point()+theme_bw()+
  ggtitle('Comparison of p-values between limmatrend and Wilcoxon test')

## ---- message = F--------------------------------------------------------
de_MAST = run_MAST(sce_clean,cl_id = "hclust_eucl",cl_ref = "B-cell",norm=T,
                   set_thresh=T,fc_cutoff = 1, alpha=0.05)

## ---- message = F--------------------------------------------------------
de_MAST = run_MAST(sce_clean,cl_id = "hclust_eucl",cl_ref = "B-cell",norm=T,
                   set_thresh=F,fc_cutoff = 1, alpha=0.5)

## ----fig.width=10--------------------------------------------------------
de_MAST = merge(de_MAST, mean_exprs_dt)

p = generic_scatterplot(de_MAST, x_col = "mean_exprs", y_col = "log2FC",
                        color = "DE_flag")
print(p+ggtitle("MA plot for MAST"))

table(de_MAST$DE_flag)
table(de_MAST[de_limma_trend$gene_id,]$DE_flag,de_limma_trend$DE_flag)

# Compare p-values between MAST and limma
ggplot(data.table(MAST = de_MAST[de_limma_trend$gene_id,on="gene_id"]$pval,
                  limmatrend = de_limma_trend$P.Val),
       aes(x=-log10(MAST),y=-log10(limmatrend)))+geom_point()+theme_bw()+
  ggtitle('Comparison of p-values between limmatrend and MAST')

## ----fig.width=10, eval = T----------------------------------------------
merged = merge(de_MAST,
               de_limma_trend[, c("gene_id","DE_flag"),with=F],
               by = "gene_id",all=T)
merged = merge(merged, de_wilcox[,c("gene_id","DE_flag"),with=F],by="gene_id",all=T)

setnames(merged,which(grepl("DE_flag",names(merged))),paste0("DE_flag_",c(1:3)))
merged[,overlap:=factor(DE_flag_1==T & DE_flag_2==T & DE_flag_3==T, labels =c("Not DE or not found by all methods","Overlap of all methods"))]

p = generic_scatterplot(merged, x_col = "mean_exprs", y_col = "log2FC", color = "overlap" )
print(p+ggtitle("Overlap between all methods tested") + ylab("log2FC (MAST)"))

## ---- eval = T, message = F----------------------------------------------
library(SC3)

# add the slots the calc biology function checks
dummy = list(`0`=0)
sce_clean@sc3$consensus = dummy
sce_clean@sc3$n_cores = 4
fData(sce_clean)$sc3_gene_filter = rep(TRUE,dim(sce_clean)[1]) #do not filter out any genes

# add whatever clustering assignment you want to the sc3_0_clusters slot
sce_clean$sc3_0_clusters = as.factor(sce_info$cell_type) 

sce_clean = sc3_calc_biology(sce_clean, k = 0)

## ---- out.width="110%",fig.height=6--------------------------------------
# change the plotted gene names to symbol for better readability
plot_sce = sce_clean[!is.na(fData(sce_clean)$symbol),]
rownames(plot_sce) = fData(plot_sce)$symbol
custom_sc3_plot_markers(plot_sce, k=0, p.val = 0.01, auroc = 0.90, show_pdata='sc3_0_clusters')
