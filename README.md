# bu2s

BU2S is a forward-time, mutation-based, infinite-sites, individual-based, stochastic population genetic simulator.
It includes mutation, selection, migration, recombination, and drift.  A full explanation of the biological scenarios simulated with this program can be found in the following two publications:

1. Flaxman SM, Wacholder AC, Feder JL, and Nosil P. 2014.  Theoretical models of the influence of genomic architecture on the dynamics of speciation.  Molecular Ecology 23:4074-4088. [doi: 10.1111/mec.12750](http://dx.doi.org/10.1111/mec.12750).
2. Flaxman SM, Feder JL, and Nosil P.  2013.  Genetic hitchhiking and the dynamic buildup of genomic divergence during speciation-with-gene-flow. Evolution 67:2577-2591. [doi: 10.1111/evo.12055](http://dx.doi.org/10.1111/evo.12055).

(see [https://www.researchgate.net/profile/Samuel_Flaxman/publications](https://www.researchgate.net/profile/Samuel_Flaxman/publications) for downloads of PDFs)

Example data sets produced by the code for the publications above are available at DataDryad.org ([http://dx.doi.org/10.5061/dryad.kc596](http://dx.doi.org/10.5061/dryad.kc596), and [http://dx.doi.org/10.5061/dryad.t894r](http://dx.doi.org/10.5061/dryad.t894r)) 

Additional capabilities since the publication of those articles includes (1) modeling neutral sites in the genome and (2) collection of additional population genetic metrics.

System Requirements for compiling and running the program:

* a Linux/UNIX-style terminal or terminal emulator
* a C compiler (gcc is the default in the Makefile)
* the [GNU Scientific Library](https://www.gnu.org/software/gsl/)
* make

To use:

(1) Download the "Source" directory from this GitHub repository.

(2) In a Linux/UNIX terminal type the following (not the dollar signs; those represent the prompt):
    
    $ cd Source
    $ make clean
    $ make

The executable generated by make ("bu2s") can then be run with the command:
    
    $ ./bu2s

A list of command line options and their explanations can be obtained with the command:
    
    $ ./bu2s -?


This work bu2s is distributed under the Creative Commons Attribution-NonCommercial 4.0 International License (CC BY-NC 4.0) [https://creativecommons.org/licenses/by-nc/4.0/](https://creativecommons.org/licenses/by-nc/4.0/).

bu2s is distributed with no warranty whatsoever.  It does not have even an implied warranty of merchantability or fitness for particular purpose.

Please email questions, comments, complaints, ideas, requests, and awesome cat pictures to samuel.flaxman@colorado.edu

