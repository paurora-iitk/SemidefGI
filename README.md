# SemidefGI

The main file in the repository is graph_iso_test.m It implements the function by the same name in the paper https://doi.org/10.1007/s10878-020-00598-w. The initial arguments are the adjacency matrices of the input graphs G_1, G_2, a n-dimensional vector, map, having all entries -1 and itr initialized to 0. The function returns a flag and the updated value of itr. The function returns flag=1 if the graphs are isomorphic, and flag=0 otherwise.

The DADAL package from https://www.math.aau.at/or/Software/ is used to solve the SDPs.

The data files are the *.mat files created with graphs taken from http://www.maths.gla.ac.uk/~es/srgraphs.php For example the srg_50_1-5.mat file contains the adjacency matrices for the first five graphs listed at http://www.maths.gla.ac.uk/~es/SRGs/50-21-8-9
