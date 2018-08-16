# NNLS
This work builds on an approach first outlined in our paper 'A Single Community Dominates Structure and Function of a Mixture of Multiple Methanogenic Communities' (see: https://doi.org/10.1016/j.cub.2017.09.056 ). The goal was to inspect the microbial community that results from mixing a number of separate (seed) communities, identify each seed community's contribution toward the final mixture and identify which contributed most to the final community. 

All the samples were sequenced (16S rRNA) and an OTU (operational taxonomic unit) counts table constructed. From this, two types of data structure were obtained: a matrix 'A' (m rows of OTUs by n sample columns) for all of the single seed communities, and a number of column vectors 'b' one for each of the mixed communities; both 'A' and 'b' hold non-negative integers of OTU abundances. The contribution, or weight, of each seed sample to the pattern of OTUs observed in each of the mixed communities is given by the column vector 'x' when solving a system of linear equations of the form 'Ax = b'.

Non-negative least squares [NNLS] to solve a system of linear equations for a non-square [m > n] matrix
