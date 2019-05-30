######################################################
# Setupfile for running SCDE on the PBMC example data
######################################################

library(scater) #needed to load the SCESets

#loading data
load("/da/dmp/cb/wegmare1/scRNASeq/data/example/output/sce_clean.RData")
load("/da/dmp/cb/wegmare1/scRNASeq/data/example/output/sce_info.RData")

# the count matrix
# Note that raw counts are required as input
counts = counts(sce_clean) 

# the cell assignments as a named factor
# here, we collapse all the non-B-cells into a single group "Others",
# because scde cannot work with very small groups
groups = sce_info$hclust_eucl
levels(groups)[groups!="B-cell"] = "Others"
names(groups) = colnames(sce_info)

# directory to save output to
out_data_dir = "/da/dmp/cb/wegmare1/scRNASeq/data/example/output/" 

# cores to use
ncores=8

# define list of groups to compare. Entries of this list are the
# identifier of the group you want to compare the rest against
diff_pair_list 		= list(
  Bcell_vs_Rest = "B-cell"
)