####################################################################
## Tutorial on how to use NNLS_experiment
##
##
#####################################################################

library(ggplot2)

# Little helper to create pooling matrix
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
# make pool
pool=makePool(3)

# read in data and create obj

myObj=read_data("CSS_norm_mcf-0.001.biom","mappingFile.txt", pool = pool)

singles=c(1,2,3)
#singles=c(1,2,3,4,5,6,7,8,9)
mixed=c(4,5,6,7,8,9,10)
#mixed=c(10,11,12,13,14,15,16,17,18,19,20,21,22,23,24,25,26,27,28,29,30)

# Set up system of lin equations

A=set_up_A(myPhyloSeq = myObj, level="Family", singles)
b=set_up_b(myPhyloSeq = myObj, level="Family", mixed)
weightSolutions=runNNLS(A, b)
write.table(weightSolutions, sep="\t", "refactored_solutionWeights.txt")
ggsave(file = "multipage_A4.pdf", outputBarcharts(weightSolutions), width = 8.27, height = 11.69, unit = "in")


