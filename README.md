# BORICE.genomic.v3
This program (Borice.Genomic.v3.py) includes a number of updates to the code in Borice.genomic (Borice.genomic.dt1.py and Borice.genomic.dt2.py).  This package should replace that one. The primary changes for version 3:

1.  The model estimates a seperate value for t (the outcrossing rate) for each population specified in the popfile.  The v3 program outputs a posterior density for each population.
2.  The v3 program allows the user to run the program without data to investigate the prior distributions implied by parameter settings specified in Control.txt
3.  The v3 program takes a replication ID as a command line argument.  This facilitates parallel use of the program (different chains run simultaneously with outputs to files specific to the replication ID).
4.  
# v3: Output posterior density for Dxy for each pairwise contrast of populations



Borice.Genomic.v3.py written in Python 2.7 (http://www.python.org/) and tested using the Linux operating system.  The program is dependent on SciPy.
If libraries are accessible, the programs should run using Python 2.7 in Windows or Mac operating systems. The command line usages:



Key parameters, such as names of input genotype file, are indentified in Control.txt..
