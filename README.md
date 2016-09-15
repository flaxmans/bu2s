# bu2s

BU2S is a forward-time, mutation-based, infinite-sites, individual-based, stochastic population genetic simulator.
It includes mutation, selection, migration, recombination, and drift.  A full explanation of the biological scenarios simulated with this program can be found in the following two publications:

1. Flaxman SM, Wacholder AC, Feder JL, and Nosil P. 2014.  Theoretical models of the influence of genomic architecture on the dynamics of speciation.  Molecular Ecology 23:4074-4088. doi: 10.1111/mec.12750.
2. Flaxman SM, Feder JL, and Nosil P.  2013.  Genetic hitchhiking and the dynamic buildup of genomic divergence during speciation-with-gene-flow. Evolution 67:2577-2591. doi: 10.1111/evo.12055.

Additional capabilities since the publication of those articles includes (1) modeling neutral sites in the genome and (2) collection of additional population genetic metrics.

System Requirements for compiling and running the program:

* a Linux/UNIX-style terminal or terminal emulator
* a C compiler (gcc is the default in the Makefile)
* make

To use:

(1) Download the "Source" directory from this GitHub repo

(2) In a Linux/UNIX terminal type the following (not the dollar signs; those represent the prompt):
    
    $ cd Source
    $ make clean
    $ make

The executable generated by make ("bu2s") can then be run with the command:
    
    $ ./bu2s

A list of command line options and their explanations can be obtained with the command:
    
    $ ./bu2s -?

