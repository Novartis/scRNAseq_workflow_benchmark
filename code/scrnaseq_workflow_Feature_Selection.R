#####################################################
#
# scRNASeq pipeline functions
#
# PART IV: Feature selection
# _______________________
#
# This script contains wrapper functions to different feature selection methods.
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


###########################################
# Feature selection
###########################################

# Based on GO terms
#_____________________________
# go_id = a GO identifier, e.g. "GO:0007049"
# returns a list of ensembl gene identifiers annotated with this GO ID or one of its child terms
GO_to_gene = function(go_id, organism = "human"){
  if(organism == "human"){
    library(org.Hs.eg.db)
    gene_eg = get(go_id,org.Hs.egGO2ALLEGS) # retrieve all entrez gene identifiers mapping to go_id
    gene_ens = unlist(lapply(gene_eg, function(x) get(x,org.Hs.egENSEMBL))) #convert to ensembl
  } else if(organism == "mouse"){
    library(org.Mm.eg.db)
    gene_eg = get(go_id,org.Mm.egGO2ALLEGS) # retrieve all entrez gene identifiers mapping to go_id
    gene_ens = unlist(lapply(gene_eg, function(x) get(x,org.Mm.egENSEMBL))) #convert to ensembl
  } else {stop("Organism has to be human or mouse.")}
  return(gene_ens)
}

# Based on highly variable genes (Brennecke 2013)
#_______________________________

info.genes = function(x,PLOT=F,qcv=.3,pv=.05,q=.95,minBiolDisp=0.1, perc_genes = NULL){
  
  if(!(is.null(perc_genes)|is.null(pv))){
    stop("Please provide either pv or perc_genes, not both!")
  }
  library(statmod)
  
  # calculate mean, variance and CV
  means  = rowMeans(x)
  vars  =  (apply(x,1,var))
  cv2 = vars/means^2
  
  # exclude genes with very low mean
  minMeanForFit  =  unname( quantile( means[ which( cv2 > qcv) ], q ) )#is the 95% of the means with a dispersion greter than 0.3
  useForFit  =  means >= minMeanForFit
  
  #fit model
  fit  =  glmgam.fit( cbind( a0 = 1, a1tilde = 1/means[useForFit] ),cv2[useForFit] ) # linear fit
  fit$coefficients
  
  a0 = unname(fit$coefficients["a0"])
  a1 = unname(fit$coefficients["a1tilde"]) #we assume xi = mean(technical size factors) = 0
  minBiolDisp = minBiolDisp^2 #squared minimum biological variance
  psia1theta = a1 # this is the term psi+a1*theta that appears in the formula for omega
  # again, we assume that the technical sf = 0 and mean ratio of all size factors = 1
  m  =  ncol(x)
  cv2th  =  a0 + minBiolDisp + a0 * minBiolDisp #a0 adjusted for min biol variation
  testDenom  =  ( means * psia1theta + means^2 * cv2th ) / ( 1 + cv2th/m ) #omega
  
  pval  =  pchisq(vars*(m-1)/testDenom,df=m-1,lower.tail=F) #Chi^2 distribution
  adj.pval  =  sort(p.adjust(pval,"fdr"))
  
  if(!is.null(pv)){
    info = adj.pval < pv
  } else {
    info = adj.pval < adj.pval[as.integer(perc_genes/100*length(adj.pval))]
  }
  
  if(PLOT){
    if(min(means)<=0) ps = .1-min(means)
    if(min(means)>0)  ps = 0
    xg  =  10^(seq( min(log(means+ps)), max(log(means+ps)), length.out=5000 ))
    
    vfit  =  a1/xg + a0 #estimated technical variation
    vfit_biol = psia1theta/xg + a0 + minBiolDisp # expected total variation
    
    xlabel = "log[ Average normalized read count]"
    smoothScatter(log(means+ps),log(cv2),xlab=xlabel,ylab="log[ Squared coefficient of variation (CV^2)]")
    points(log(means+ps),log(cv2),col="gray")
    # lines(log(xg[which(vfit>0)]+ps), log(vfit[which(vfit>0)]), col="black", lwd=3,lty=2,ylim=range(cv2) )
    lines(log(xg[which(vfit_biol>0)]+ps), log(vfit_biol[which(vfit_biol>0)]), col="black", lwd=3,ylim=range(cv2) )
    lines(log(xg[which(vfit_biol>0)]+ps),log(vfit_biol[which(vfit_biol>0)] * qchisq(0.975,m-1)/(m-1)),lty=2,col="black")
    lines(log(xg[which(vfit_biol>0)]+ps),log(vfit_biol[which(vfit_biol>0)] * qchisq(0.025,m-1)/(m-1)),lty=2,col="black")
    points(log(means[names(which(info)[TRUE])]+ps),log(cv2[names(which(info)[TRUE])]),col=2,cex=.5)
    
  }
  
  return(info)
}


# M3Drop, using the two DANB based methods
#_______________________________________________

run_DANB = function(counts,save_plot=T,method="NBDrop",cutoff=NULL, perc_genes = NULL){
  
  if((is.null(perc_genes)&is.null(cutoff))){
    stop("Please provide exactly one of either cutoff or perc_genes")
  }
  
  if(!(is.null(perc_genes)|is.null(cutoff))){
    stop("Please provide either cutoff or perc_genes, not both!")
  }
  
  library(M3Drop) #note that you require version > 2.0 (not on bioconductor yet)
  
  if(!method%in%c("NBDrop","NBDisp")){
    stop("Invalid method selected. Please choose one of \"NBDrop\", \"NBDisp\"")
  }
  
  # fitting the dropout model (DANB)
  fit = NBumiFitModel(counts) #this fits a DANB model
  fit$mus = t(sapply(fit$vals$tjs, function (x) x * fit$vals$tis/fit$vals$total))
  
  size_coeffs = NBumiFitDispVsMean(fit, suppress.plot = T)#get coefcients of mean-dispersion fit
  smoothed_size = exp(size_coeffs[1] + size_coeffs[2] * log(fit$vals$tjs/fit$vals$nc)) #predicted dispersions per gene
  size_mat = matrix(rep(smoothed_size, times = fit$vals$nc), ncol = fit$vals$nc, byrow = FALSE)
  exp_ps  =  (1 + fit$mus/size_mat)^(-size_mat) #these are the fitted values per cell and gene
  exp_tot = rowSums(exp_ps) #this is the total predicted molecules per gene
  
  plot_dt = data.table(Dropout_rate = fit$vals$djs/fit$vals$nc,
                       expression = fit$vals$tjs/fit$vals$nc,
                       sizes = fit$sizes,
                       pred_sizes = smoothed_size,
                       predicted_dropouts=exp_tot/fit$vals$nc)
  
  if(method=="NBDrop"){
    NBumiCheckFitFS(counts,fit,suppress.plot = T) #check whether the fitted droout rates are well correlated with observed ones (i.e. number of zeroes)
    pvals = NBumiFeatureSelectionCombinedDrop(fit) #ranks genes by difference from expected dropout rates
    
    if(is.null(cutoff)){ cutoff = sort(pvals)[as.integer(perc_genes/100*length(pvals))]}
    
    info = rownames(counts)%in%names(pvals[which(pvals<cutoff)]) #select the top features based on dropout rates
    plot_dt[,info:=info]
    
    p = ggplot(plot_dt,aes(x=log10(expression),y=Dropout_rate)) +
      geom_point(aes(color=info),alpha=0.7,size=2) + 
      geom_line(aes(x=log10(expression),y=predicted_dropouts),color=colors()[30],size=1.2)+
      theme_bw() + xlab("Expression [log10]") + ylab("Dropout Rate")+
      ggtitle("Dropout rate vs. Expression")+ theme(text = element_text(size=17))+
      scale_color_manual(values = colors()[c(226,32)],name="is_outlier")
    print(p)
    if(save_plot){
      ggsave(p,file=file.path(plotdir,"NBdrop_plot.pdf"),height=6,width=8)
    }
  } else {
    resids = NBumiFeatureSelectionHighVar(fit) #ranks genes by difference from fitted mean-dispersion power law
    if(is.null(cutoff)){cutoff = perc_genes/100}
    info = rownames(counts)%in%names(resids[which(resids<quantile(resids,cutoff))])
    plot_dt[,info:=info]
    plot_dt[,est_cv2:=(1+expression/sizes)/expression] #predicted variance according to DANB model
    plot_dt[,ave_cv2:=(1+expression/pred_sizes)/expression] #predicted variance according to linear fit of dispersions
    
    p = ggplot(plot_dt,aes(x=log10(expression),y=log10(est_cv2))) +
      geom_point(aes(color=info),alpha=0.7,size=2) + 
      geom_line(aes(x=log10(expression),y=log10(ave_cv2)),color=colors()[30],size=1.2)+
      theme_bw() + xlab("Expression [log10]") + ylab("Estimated CV^2 [log10]")+
      ggtitle("Mean - Dispersion Relationship")+ theme(text = element_text(size=17))+
      scale_color_manual(values = colors()[c(226,32)],name="is_outlier")
    print(p)
    if(save_plot){
      ggsave(p,file=file.path(plotdir,"NBdisp_plot.pdf"),height=6,width=8)
    }
  }
  return(info)
}
