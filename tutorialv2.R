####################################################################
## Tutorial on how to use NNLS_experiment
##
##                                              
#####################################################################

library(ggplot2)

# make pool
pool=makePool(3)

# read in data and create obj

myObj=read_data("NNLS/CSS_norm_mcf-0.001.biom","NNLS/mappingFile.txt", pool=pool)

singles=c(1,2,3)
mixed=c(4,5,6,7,8,9,10)

# Set up system of lin equations

A=set_up_A(myPhyloSeq = myObj, level="Family", singles)
b=set_up_b(myPhyloSeq = myObj, level="Family", mixed)

write.table(runNNLS(A, b), sep="\t", "refactored_solutionWeights.txt")
weightSolutions=runNNLS(A, b)
ggsave(file = "multipage_A4.pdf", outputBarcharts(weightSolutions,myObj), width = 8.27, height = 11.69, unit = "in")


