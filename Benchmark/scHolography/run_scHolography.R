### Installation
if (!requireNamespace("remotes", quietly = TRUE)) {
  install.packages("remotes")
}
remotes::install_github(sprintf("rstudio/%s", c("reticulate", "tensorflow", "keras")))
reticulate::miniconda_uninstall() 
reticulate::install_miniconda()
tensorflow::install_tensorflow(method = 'conda',envname = 'r-tensorflow',version = '2.12',conda_python_version = '3.9')

######### After above steps, for Mac Apple silicon users, you need to install tensorflow-macos instead of tensorflow by first activating the conda environment and then installing tensorflow-macos in your terminal
## Used command 'conda env list' to find the environment path, like /home/yuling/.local/share/r-miniconda/envs/r-tensorflow
## conda activate /home/yuling/.local/share/r-miniconda/envs/r-tensorflow
## python -m pip install tensorflow-macos==2.12.0
######### Then you can set the RETICULATE_PYTHON to the python in the conda environment in your R/RStudio
## Sys.setenv(RETICULATE_PYTHON = "/Users/YOURPATH/Library/r-miniconda/envs/r-tensorflow/bin/python")

reticulate::use_condaenv('r-tensorflow')
library(tensorflow)
tf$constant("Hello Tensorflow!")
options(repos = c("https://satijalab.r-universe.dev", getOption("repos")))
remotes::install_version("SeuratObject", "4.1.3")
remotes::install_version("Seurat", "4.3.0", upgrade = FALSE) 
remotes::install_github("YiLab-SC/scHolography", upgrade = FALSE)
### Load Package
library(scHolography)
library(Seurat)
library(ggplot2)
library(dplyr)
library(RColorBrewer)
library(viridis)
library(Gmisc)
library(dplyr)
library(plotly)
## Load mHypothalamus data both for reference data (low.res.sp) and query data (high.res.sp) 
create_obj<-function(data_name,slice_order){
  matrix_dir=pathJoin('/home/data3/yuling/xenium_region2_5599/Dataset/matrix',data_name)
  mHypothalamus_9<-Read10X(data.dir = matrix_dir)
  metadata<-read.csv(pathJoin('/home/data3/yuling/xenium_region2_5599/Dataset/matrix',data_name,'metadata.csv'))
  rownames(metadata)<-metadata$X
  sobj<-CreateSeuratObject(counts=mHypothalamus_9,meta.data = metadata,assay = 'RNA')
  sobj@meta.data$slice=data_name
  return (sobj)
}
#### reference data needs to have spatial information
data_list=c('-0.04','-0.09','-0.14','-0.19','-0.24')
data_1<-create_obj(data_name = data_list[1],slice_order = 1)
data_1<-subset(data_1,cells=which(is.na(data_1$original_clusters)==F))
data_1@meta.data[,c("CentroidX","CentroidY")]<-data_1@meta.data[,c("x","y")]
com_obj<-data_1
data_1@images$image =  new(
  Class = 'SlideSeq',
  assay = "Spatial",
  key = "image_",
  coordinates = data_1@meta.data[,c("CentroidX","CentroidY")]
)
for (i in 2:5){
  data_2<-create_obj(data_name = data_list[i],slice_order = i)
  data_2<-subset(data_2,cells=which(is.na(data_2$original_clusters)==F))
  # align data to reference slice
  data_2@meta.data[,c("CentroidX","CentroidY")]<- t(t(data_2@meta.data[,c("x","y")])-c(median(data_2@meta.data[,c("x")]),median(data_2@meta.data[,c("y")]))+c(median(data_1@meta.data[,c("x")]),median(data_1@meta.data[,c("y")])))

    com_obj<-merge(
      x = data_2,
      y = com_obj,
      project = "CombinedProject")
}
# if high.res.sp data has not been preprocessed, set scProcessed = F
# if low.res.sp data has not been preprocessed, set stProcessed = F
sp.integrated <- dataAlign(low.res.sp =data_1,high.res.sp =com_obj,nPCtoUse = 32,scProcessed = F,future.size = 8000)
scHolography.obj<-trainHolography(sp.integrated,n.slot = 30,n.pcUse = 32,n.pcOut = 32,n.repeat = 30)

### Set the original cluster as factor
scHolography.obj$scHolography.sc$original_clusters<- factor(as.character(scHolography.obj$scHolography.sc$original_clusters),levels = c('MPA','MPN','BST','fx','PVH','PVT','V3','PV'))

### To orient the plot in a more intuitive way
scHolography.obj$scHolography.sc$z3d_sp <- -scHolography.obj$scHolography.sc$z3d_sp
scHolography.obj$scHolography.sc$x3d_sp <- -scHolography.obj$scHolography.sc$x3d_sp
scHolography.obj$scHolography.sc$y3d_sp <- -scHolography.obj$scHolography.sc$y3d_sp
scene = list(camera = list(eye = list(x = -1, z = 0, y = 2)))

### To visualize the 3D structure colored by celltype
scHolographyPlot(scHolography.obj,color.by = "original_clusters")%>% plotly::layout(scene = scene)
### cluster distance boxplot
clusterDistanceBoxplot(scHolography.obj,annotationToUse = "original_clusters",reference.cluster = "MPA",query.cluster.list = c('BST','MPA','MPN','PV','PVH','PVT','V3','fx'))
### scHolographyNeighborComPlot 
neighbor.comp<- scHolographyNeighborCompPlot(scHolography.obj,annotationToUse = "original_clusters",query.cluster = c('BST','MPA','MPN','PV','PVH','PVT','V3','fx'))

my.color.order= c('#1F77B4', '#FF7F0E', '#2CA02C', '#D62728', '#9467BD', '#8C564B', '#E377C2', '#7F7F7F') # color order for the plot
neighbor.comp$neighbor.comp.plot+scale_fill_manual(values = my.color.order)
neighbor.comp$neighbor.comp.sig.plot+scale_fill_manual(values = my.color.order)
