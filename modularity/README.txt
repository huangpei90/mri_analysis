This is the Matlab code used for calculating network modularity, recruitment and intergration using the community louvain function.
Requires BCT toolbox (download here: https://sites.google.com/site/bctnet/)

Input: 
1. Subject directory (in_dir) with connectivity matrices arranged in running subject order
2. comb_in and comb_out according to the grouping of ROIs into networks (in example, ROIs 1-6 belong to RP network, 7-10 belong to CC network, and 11-14 belong to EP network)

Output: 
1. Tabulated results (outfile)
