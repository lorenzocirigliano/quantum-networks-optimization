This folder contains the program "average_optimal_networks", which can be used to do numerical simulations and compute the average properties of the optimal network and reproduce the results in the paper.

The program requires as input:

1) N: number of nodes
2) L/lambda_0: rescaled length of the square, 
3) p: probability of node leakage, 
4) alpha_MIN: minimum value of alpha for which you want to optimize, 
5) alpha_MAX: maximum value of alpha for which you want to optimize, 
6) Delta_alpha: the resolution you want for alpha, 
7) M: the number of realizations on which you want to average,
8)BOOLE: a binary {0,1} variable: set BOOLE = 1 if you want to compute degree distribution, link density and betweenness. This slows down the execution of the program. Set BOOLE = 0 if you don't want to compute the average topological properties of the optimal networks. 

After having compiled the .c file, make executable the file average_optimal_networks.sh (in this file you can change the parameters as you prefer).
Then it is sufficient to run the command

./average_optimal_networks.sh
