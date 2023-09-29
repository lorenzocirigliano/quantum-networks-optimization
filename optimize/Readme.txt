This folder contains the program that, given the coordinates of N trusted nodes,
produces in output the list of connections of the network with optimal efficiency.

Please proceed as follows.



*** (1) ***
Compile the file optimize_QE.c
producing an executable optimize_QE

For example, for Linux users the command line is:

gcc optimize_QE.c -lm -o optimize_QE



*** (2) ***
Modify the script optimize.sh (make it executable if needed)
to specify the parameters needed in input:

- The name <input_file_name> of a file containing two columns with the
  spatial coordinates (x,y) of the N trusted nodes.
- The values of the following parameters:
  N = number of nodes
  lambda_0 = decay length of the capacitance
  p = probability that a node is not secure
  Delta_alpha = interval between subsequent values of alpha considered



*** (3) ***
Run the script

./optimize.sh

and follow the instructions. The procedure is divided in two parts:

- In Part 1 the script performs the optimization algorithm for values
  of alpha between 0 and 1 in steps Delta_alpha and
  produces an output file "optimal_efficiency_for_<input_file_name>.dat"
  with three columns, containing:
  alpha   Q*  L* 

  The program stops and waits for a value of alpha.
  Please inspect the output file, decide which value of (Q*, L*) you
  want and check the corresponding value of alpha.
  Insert this value of alpha at the prompt.

- In Part 2 the program creates the optimal network using the modified Pollack algorithm.
  The output is the file "list_of_connections_N=*_lambda_0=*_p=*_alpha=*_for_<input_file_name>.dat",
  containing in each row the labels of the nodes that must be connected by a link.
  The node labels are ordered from 0 to N-1 according to their ordering in <input_file_name>.

