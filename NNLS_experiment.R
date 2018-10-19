####################################################################
## Solving a system of linear equations for a non-square          ##
## matrix via the use of non-negative least squares [NNLS].       ##
## Returns the solution 'weights' and residual sum of squares.    ##
##                                                                ##
## Mark Alston,                                                   ##
## Cell & Developmental Biology,                                  ##
## John Innes Centre, Norwich, UK.
## mark.alston@jic.ac.uk
##
## Sarah Bastkowski,
## Earlham Institute, Norwich, UK.
## sarah.bastkowski@earlham.ac.uk
##
##
#####################################################################



## ========================================================================================================= ##
## READ IN THE OTU table
## ========================================================================================================= ##
## MAPPING FILE TO BE USED: 'mappingFile.txt'
## OTU TABLE TO BE USED:    'CSS_norm_rawValues.biom'



#Set up function that reads in files, pooles samples

read_data <- function(biomfile, mappingfile, pool=NULL, cellcounts=NULL) {
  
  library("biomformat")
  library("phyloseq")
  
  data = import_biom(biomfile)
  SAMPLES = import_qiime_sample_data(mappingfile)
  myTaxTable <- tax_table(data)
  colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
  OTU = otu_table(data)
  TAX = tax_table(myTaxTable)
  myPhyloSeq <- phyloseq(OTU,TAX)
  
  if(!missing(cellcounts)){
    myPhyloSeq=normByCellcount(myPhyloSeq, cellcountvalues = cellcounts)
  }
  
  if(!missing(pool)){
    OTU_pooled=poolSamples(myPhyloSeq, pool)
    myPhyloSeq <- merge_phyloseq(OTU_pooled, TAX)
  }
  
 
  
  return(myPhyloSeq)
}

#Normalise by cellcount

normByCellcount <- function(phyloseqObj, cellcountvalues) {
  pooledOTU=rep(0,nrow(otu_table(phyloseqObj)))
  
  for(i in 1:length(cellcountvalues)){
    for(j in 1:nrow(phyloseqObj@otu_table@.Data)){
      phyloseqObj@otu_table@.Data[j,i]=phyloseqObj@otu_table@.Data[j,i]*cellcountvalues[i]
    }
  }
  
  return(phyloseqObj)
}

# Input: OTU table and matrix that indicates which samples should be pooled,
# row is one pool and columns (number of samples) indicates which Samples are in this pool
# if entry is 0 its not in pool, if sample is 1 it is in pool.
# single samples are treated as a pool as well, ie one row with zeros except one 1
poolSamples<-function(phyloseqObj, pool) {
  pooledOTU=rep(0,nrow(otu_table(phyloseqObj)))
  
  for(i in 1:nrow(pool)){
    C=matrix(rep(0, nrow(otu_table(phyloseqObj))))
    for(j in 1:ncol(pool)) {
      C=C+ pool[i,j]*otu_table(phyloseqObj)[,j]
    }
    pooledOTU=cbind.data.frame(pooledOTU,C/sum(pool[i,]))
    
  }
  row.names(pooledOTU)=row.names(tax_table(phyloseqObj))
  pooledOTU=pooledOTU[,-1]
  
  return(otu_table(pooledOTU, taxa_are_rows = TRUE))
}


## ========================================================================================================= ##
## collapse OTU table at <level> given as string, eg "Family", and generate matrix 'A' (single communities)
## and 'b' (mixed communities)
## Input: the phyloseq object, level, vector stating which columns in the matrix (OTU) are single samples (A)
## and which ones are mixed (b)
## ========================================================================================================= ##
set_up_A <- function(myPhyloSeq, level, singles, mixed) {
  #Combine A and b in dataframe including samplenames
  bacteria_level <- tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(get_taxa(otu_table(bacteria_level)))
  #need to extract this
  bacteria_level_singleComm_df  <-  bacteria_level_df[,singles]
  # A
  bacteria_family_matrix_A <- data.frame(as.matrix(bacteria_level_singleComm_df))
  colnames(bacteria_family_matrix_A)=colnames(myPhyloSeq@otu_table@.Data)[singles]
  
  return(bacteria_family_matrix_A)
}

set_up_b <- function(myPhyloSeq, level, mixed) {
  # b
  my_bvectors <- list()
  bacteria_level <- tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(get_taxa(otu_table(bacteria_level)))
  for (i in 1:length(mixed))
  {
    my_bvectors[[i]] <-bacteria_level_df[,mixed[i]]
  }
  
  names(my_bvectors)=colnames(myPhyloSeq@otu_table@.Data)[mixed]
  return(my_bvectors)
  
}

runNNLS<-function(matrix, vector){
  library("nnls")
  my_weightsList <- list()
  
  for (i in 1:length(vector))
  {
   
    weightObject <-  nnls(as.matrix(matrix),vector[[i]])
    
    my_weightsList[[i]] <- data.frame(weightObject$x, colnames=colnames(matrix))
    names(my_weightsList)[i]=names(vector)[i]
    }
  #add names for samples to be used in barchart
 
  return(my_weightsList)
}


weighttoPercent<-function(weights){
  
  percentages=list()
  for (i in 1:length(weights))
  {
    percent=vector(mode="integer", length = length(weights[[i]]$weightObject.x))
    for (j in 1:length(weights[[i]]$weightObject.x)){
      percent[j] <-  weights[[i]]$weightObject.x[j]/sum(weights[[i]]$weightObject.x)
    }
    percentages[[i]]=percent
  }
  #add names for samples to be used in barchart
  
  return(percentages)
}
## ========================================================================================================= ##
## THIRD: Plot, arrange and output the solution 'weight' values as barcharts
## ========================================================================================================= ##


######## THE BARCHART OUTPUT FUNCTION ##############

outputBarcharts <- function(weightSolutions, percent=TRUE){
  library(lattice)
  library(gridExtra)
  library(ggplot2)
  
  if (percent==TRUE){
    p=weighttoPercent(weights = weightSolutions)
    trellis.objects.list = list()
    
    for (i in 1:length(weightSolutions))
    {
      trellis.objects.list[[i]] <- barchart(as.table(p[[i]]),
                                            main=names(weightSolutions)[i],
                                            horizontal=FALSE, col="steelblue", ylab="Weight", xlab = "Samples", aspect=1,
                                            scales=list(x=list(rot=70, labels=weightSolutions[[i]]$colnames, cex=1.1)))
    }
    
    
    # Output and arrange the generated barchart
    # see: https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
    # Need to create graphical objects (grobs) to use 'grid.arrange' etc. ...whatever...
    
    # multi-page output probably best for most users
    multiPageGrobs <- marrangeGrob(grobs = trellis.objects.list, nrow=4, ncol=2)
    
  }else{
    # create a list to hold the barchart "trellis" objects
    trellis.objects.list = list()
    
    for (i in 1:length(weightSolutions))
    {
      trellis.objects.list[[i]] <- barchart(as.table(weightSolutions[[i]]$weightObject.x),
                                            main=names(weightSolutions)[i],
                                            horizontal=FALSE, col="steelblue", ylab="Weight", xlab = "Samples", aspect=1,
                                            scales=list(x=list(rot=70, labels=weightSolutions[[i]]$colnames, cex=1.1)))
    }
    
    
    # Output and arrange the generated barchart
    # see: https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
    # Need to create graphical objects (grobs) to use 'grid.arrange' etc. ...whatever...
    
    # multi-page output probably best for most users
    multiPageGrobs <- marrangeGrob(grobs = trellis.objects.list, nrow=4, ncol=2)
    
  }
 
  # Send to A4-sized output
  # note: 'ggsave' recognises the extensions eps/ps, tex (pictex), pdf, jpeg, tiff, png, bmp, svg and wmf (windows only).
  
  return(multiPageGrobs)
}
# end of 'outputBarcharts' function
