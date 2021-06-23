This is the Readme file for the Omni-Optimizer source code.


About the Algorithm
--------------------------------------------------------------------------
Omni-Optimizer

Please refer to the following paper for details about the algorithm:

Authors: Kalyanmoy Deb and Santosh Tiwari
Paper Title: Omni-Optimizer: A Generic Evolutionary Algorithm for Single and Multi-objective Optimization
Journal: European Journal of Operational Research
Year: 2008
Volume: 185
Number: 3
Month: March
Pages: 1062-1087
---------------------------------------------------------------------------


---------------------------------------------------------------------------
NOTE: This archive contains routines for ploting the variables and objectives using gnuplot. The code has been written for posix compliant operating systems and uses standard piping method provided by the GNU C library. The routines should work on any unix and unix like OS having gnuplot installed and which are posix compliant.
---------------------------------------------------------------------------


How to compile and run the program
---------------------------------------------------------------------------
Makefile has been provided for compiling the program on linux (and unix-like) systems. Edit the Makefile to suit your need. By default, provided Makefile attempts to compile and link all the existing source files into one single executable.

Name of the executable produced is: omniopt

To run the program, type: ./omniopt random_seed

Here random_seed is a real number in (0, 1) which is used as a seed for random number generator.

You can also store all the input data in a text file and use a redirection operator to give the inputs to the program in a convenient way.

You may use the following syntax: ./omniopt random_seed <inp_file.in, where "inp_file.in" is the file that stores all the input parameters.
---------------------------------------------------------------------------


About the output files
---------------------------------------------------------------------------
initial_pop.out: This file contains the initial population.

final_pop.out: This file contains the final population.

nondominated_pop.out: This file contains the final feasible population. This is the optimized result.

output.out: This file containts the data of populations at all generations.

problem_data.out: This file contains the information about input parameters read by the program.

Other output files produced by the algorithm may be discarded.
---------------------------------------------------------------------------


About the input parameters
---------------------------------------------------------------------------
See the template input file: template.in

There are many example input files provided in the sampleinput directory.
---------------------------------------------------------------------------


Defining the Test Problem
---------------------------------------------------------------------------
Check the source files sampleprob1.c, sampleprob2.c, ..., sampleprob5.c about how to define your test problem. 76 different test problems are provided as examples on how to define your own objective and constraint functions. You can also link other source files with the code depending on your need.

Check the file pdef.h: It defines the test problem that will be used by the omni-optimizer. The "#define xxxxxx" line tells the test problem to be used. If you add your own test problem, code it inside #ifdef---#endif directives, and define the corresponding string in pdef.h

The following points are to be kept in mind while writing objective and constraint functions.

1.
The code has been written for minimization of objectives (min f_i). If you want to maximize a function, you may use negetive of the function value as the objective value.

2.
A solution is said to be feasible if it does not violate any of the constraints. Constraint functions should evaluate to a quantity greater than or equal to zero (g_j >= 0), if the solution has to be feasible. A negative value of a constraint means, it is being violated.

3.
If there are more than one constraints, it is advisable (though not mandatory) to normalize the constraint values by either reformulating them or dividing them by a positive non-zero constant, so that all the constraints have values in a similar range.

It is strongly recommended that you use pdef.c to define your own test problem.

Also, if you wish to do some preprocessing (before optimization starts), write that code in pdefinit.c inside the pdefinit function. This function is called before the OmniOptimizer algorithm starts. You can use this function to load initial data from a file, perform some pre-processing etc.
---------------------------------------------------------------------------


About the files
---------------------------------------------------------------------------
global.h: Header file containing declaration of global variables and functions

rand.h: Header file containing declaration of variables and functions for random number generator

pdef.h: Header file containing the definition of the test problem that will be used by the optimizer.

allocate.c: Memory allocation and deallocation routines

crossover.c: Routines for real and binary crossover

crowddist.c: Crowding distance assignment routines

decode.c: Routines to decode binary variables

display.c: Routines to display the data in real-time using gnuplot

dominance.c: Routines to perform non-domination checking

eval.c: Routines to evaluate constraint violation

fillnds.c: Non-dominated sorting based selection

ga.c: Implementation of main function and the OmniOptimizer algorithm

global.c: Declaration of all the global variables

initialize.c: Routine to perform random initialization of population members

input.c: Routine to read the inputs provided by the user

list.c: A custom doubly linked list implementation

merge.c: Routine to merge two population into one larger population

mutation.c: Routines for real and binary mutation

output.c: Routine to print the simulation information in the output files

pdef.c: Source file that can be used to add user-provided test problems

pdefinit.c: User-provided pre-processing code

rand.c: randome number generator routines

rank.c: Rank assignment routines

report.c: Routines to write the population information in different files

sampleprob1.c: Sample multi-objective test problems

sampleprob2.c: Sample single-objective test problems

sampleprob3.c: Sample constrained single-objective test problems

sampleprob4.c: Sample problems from non-standard sources

sampleprob5.c: Few more single-objective test problems

sort.c: Randomized quick sort implementation

tourselect.c: Tournament selection routine

utility.c: auxiliary routines (not part of the algorithm)
---------------------------------------------------------------------------

Please feel free to send questions/comments/doubts/suggestions/bugs etc. to deb@iitk.ac.in
---------------------------------------------------------------------------
