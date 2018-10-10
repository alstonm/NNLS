#####################################################################
## Solving a system of linear equations for a non-square [m > n]   ##
## matrix via the use of non-negative least squares [NNLS].        ##
## Returns the solution 'weights' and residual sum of squares.     ##
##                                                                 ## 
## Mark Alston,                                                    ##
## Cell & Developmental Biology,                                   ##
## John Innes Centre, Norwich, UK.     
## mark.alston@jic.ac.uk  
## 
## Sarah Bastkowski
## Earlham Institute
## sarah.bastkowski@earlham.ac.uk
##
##                                              
#####################################################################


## ========================================================================================================= ##
## Install the R packages phyloseq, biomformat, nnls and limSolve.
## Now load as required...  
## ========================================================================================================= ##




## ========================================================================================================= ##
## READ IN THE OTU table
## ========================================================================================================= ##
## MAPPING FILE TO BE USED: 'mappingFile.txt' 
## OTU TABLE TO BE USED:    'CSS_norm_rawValues.biom' 



#Set up function that reads in files, pooles samples

read_data <- function(biomfile, mappingfile, pool) {
  
  library("biomformat")
  library("phyloseq")
  
  data = import_biom(biomfile)   
  SAMPLES = import_qiime_sample_data(mappingfile)
  myTaxTable <- tax_table(data)
  colnames(myTaxTable) <- c("Kingdom","Phylum", "Class", "Order", "Family", "Genus","Species")
  OTU = otu_table(data)
  TAX = tax_table(myTaxTable)
  myPhyloSeq_allData <- phyloseq(OTU,TAX)  
  OTU_pooled=poolSamples(myPhyloSeq_allData, pool)
  myPhyloSeq <- merge_phyloseq(OTU_pooled, TAX)
  return(myPhyloSeq)
}





# Change to samplenames 
# Need to specify how to present the pooling or add function to pull this out of mapping file
makePool<-function(size){
  
  pool=matrix(rep(0,300), ncol=30, nrow=10)
  for(i in 1:nrow(pool)){
    for(j in 1:size){
      pool[i,((i-1)*size + j)]=1
    }
  }
  return(pool)
}


#Input: OTU table and matrix that indicates which samples should be pooled, 
#row is one pool and columns (number of samples) indicates which Samples are in this pool
# if entry is 0 its not in pool, if sample is 1 it is in pool.
#single samples are treated as a pool as well, ie one row with zeros except one 1
poolSamples<-function(phyloseqObj, pool) {
  pooledOTU=rep(0,nrow(otu_table(phyloseqObj)))
 
  for(i in 1:nrow(pool)){
    C=matrix(rep(0, nrow(otu_table(phyloseqObj))))
    for(j in 1:ncol(pool)) {
      C=C+ pool[i,j]*otu_table(phyloseqObj)[,j]
    }
    pooledOTU=cbind.data.frame(pooledOTU,C)
    
  }
  row.names(pooledOTU)=row.names(tax_table(phyloseqObj))
  pooledOTU=pooledOTU[,-1]
  
  return(otu_table(pooledOTU, taxa_are_rows = TRUE))
}
singles=c(1,2,3)
mixed=c(4,5,6,7,8,9,10)
## ========================================================================================================= ##
## collapse OTU table at level given as string, eg "Family", and generate matrix 'A' (signle communities)
## and 'b' (mixed communities)
## ========================================================================================================= ##
set_up_system_of_equations <- function(myPhyloseq, level, singles, mixed) {
  
  bacteria_level <- tax_glom(myPhyloSeq, taxrank=level)
  bacteria_level_df <- as.data.frame(get_taxa(otu_table(bacteria_level)))
  #need to extract this
  bacteria_level_singleComm_df  <-  bacteria_level_df[,singles] 
  # A
  bacteria_family_matrix_A <- as.matrix(bacteria_level_singleComm_df)
  
  # b
  my_bvectors <- list()
  for (i in 1:length(mixed)) 
  {
    # create a name...
    name <- paste("b_M0",i,"_bac_family_df", sep = "")
    
    # ...and assign something to an object with that name
    bvectorObject <- assign(name, bacteria_level_df[,mixed[i]])
    
    # create a name for the b-vector
    bvector_name <- paste("bv_M0",i,sep = "")
    
    # assign something to an object with that name & add the b-vector to a list 'my_bvectors'
    my_bvectors[[i]] <- assign(bvector_name, as.matrix(bvectorObject))
    
    
  }
}




row.names(bacteria_family_singleComm_df) <- NULL
colnames(bacteria_family_singleComm_df) <- NULL


## ========================================================================================================= ##
## get matrix 'A' 
## ========================================================================================================= ##





## define which otu_table columns correspond to the MIXED communities 

#mixedSamples <- c(11,12,10,19,15,16,13,14,8,9)
mixedSamples <- c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)
mixedPooledSamples<-c(4,5,6,7,8,9,10)
#mixedSamples <- c(10,11,12)
mixedSamples
no_of_mixedSamples <- length(mixedPooledSamples)




runNNLS<-function(matrix, vector){
  library("nnls")
  my_weightsList <- list()
  
  for (i in 1:length(my_bvectors)) 
  {
    # create a name...
    weightName <- paste("weight_M",i, sep = "")
    weightObject <-  nnls(bacteria_family_matrix_A, my_bvectors[[i]])
    my_weightsList[[i]] <- assign(weightName, weightObject )
  }
  
  
  return(my_weightsList)
}


write.table(runNNLS(bacteria_family_matrix_A, my_bvectors), sep="\t", "refactored_solutionWeights.txt")
weightSolutions=runNNLS(bacteria_family_matrix_A, my_bvectors)


###### ALL OK UP TO THIS POINT ###########


# need to set the names within the lists
## name columns and rows ##
###########################



## ========================================================================================================= ##
## THIRD: Plot, arrange and output the solution 'weight' values as barcharts
## ========================================================================================================= ##






######## THE BARCHART OUTPUT FUNCTION ##############

outputBarcharts <- function(weightSolutions){
  library(lattice)
  library(gridExtra)
  library(ggplot2)
  
  # create a list to hold the barchart "trellis" objects
  trellis.objects.list = list()
  seed_sample_names=c("1","2","3")
  
  for (i in 1:length(weightSolutions)) 
  {
    trellis.objects.list[[i]] <- barchart(as.table(weightSolutions[[i]]$x),
                                          main=colnames(weightSolutions)[i],
                                          horizontal=FALSE, col="steelblue", ylab="Weight", aspect=1,
                                          scales=list(x=list(rot=70, labels=seed_sample_names, cex=1.1)))
  }
  
  
  # Output and arrange the generated barchart
  # see: https://cran.r-project.org/web/packages/gridExtra/vignettes/arrangeGrob.html
  # Need to create graphical objects (grobs) to use 'grid.arrange' etc. ...whatever...
  
  # multi-page output probably best for most users
  multiPageGrobs <- marrangeGrob(grobs = trellis.objects.list, nrow=4, ncol=2)
  
  # Send to A4-sized output
  # note: 'ggsave' recognises the extensions eps/ps, tex (pictex), pdf, jpeg, tiff, png, bmp, svg and wmf (windows only).
  
  
  return(multiPageGrobs)
} 
# end of 'outputBarcharts' function
# end of 'outputBarcharts' function
ggsave(file = "multipage_A4.pdf", outputBarcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")
