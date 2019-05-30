###############################################
# scRNASeq workflow - Functions
#
# This script sources all parts of the scRNASeq workflow.
###############################################

source(file.path(code_dir, "scrnaseq_workflow_Setup.R"))
source(file.path(code_dir, "scrnaseq_workflow_QC.R"))
source(file.path(code_dir, "scrnaseq_workflow_Normalization.R"))
source(file.path(code_dir, "scrnaseq_workflow_Feature_Selection.R"))
source(file.path(code_dir, "scrnaseq_workflow_Clustering.R"))
source(file.path(code_dir, "scrnaseq_workflow_CellSIUS.R"))
source(file.path(code_dir, "scrnaseq_workflow_DE.R"))
source(file.path(code_dir, "scrnaseq_workflow_Plotting.R"))
