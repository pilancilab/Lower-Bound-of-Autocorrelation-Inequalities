# Minimum-Overlap
Code for optimization programs and feasible solutions.

The author of this repository recommends first familiarizing oneself with the basic ideas in https://arxiv.org/abs/2201.05704: "Erdos' minimum overlap problem".

In this repository you will find code for four programs, all written in OPL and meant to be run on IBM's CPLEX software. A user may be able to somewhat easily convert between this language and another of their choosing. The programs are:

### MO_simple_primal.txt

This is the program in Section 4 of Erdos' minimum overlap problem. 

### MO_simple_dual.txt

This is the dual program of the above. By adjusting inputs and constraints, one can find strictly feasible solutions to find fairly strong bounds on the minimum overlap problem. This program will print the feasible solution to a separate txt file, along with the values of (slack in) constraints. 

### MO_convex_primal.txt

This is the program in Section 5 of Erdos' minimum overlap problem. The inputs and structure exactly match the program described. 

### MO_convex_dual.txt

This is the program is the dual program to the above, when formulated as a SOCP. Details of the program are in the appendices of Erdos' minimum overlap problem. This program will print the feasible solution to a separate txt file, along with the values of (slack in) constraints. 

## Data

Feasible solutions generate by MO_convex_dual.txt are in the folder Feasible_Solutions. The solutions given match the ones referenced in Erdos' minimum overlap problem. 

### Table 2: Dual program output I

The data of this table is in folders 10000_4000_10 and 20000_5000_10. Folder 10000_4000_10 has data corresponding to the first 15 rows of this table where (N,T,R) = (10000,4000,10). Folder 20000_5000_10 has data corresponding to the final 3 rows of this table where (N,T,R) = (20000,5000,10).


### Table 3: Dual program output II

The data of this table is in folders 25000_7000_10. The names of the files are in accordance with the colours described in the table. 


