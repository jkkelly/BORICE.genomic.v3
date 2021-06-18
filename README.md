# BORICE.genomic.v3
This program (Borice.Genomic.v3.py) includes a number of updates to the code in Borice.genomic (Borice.genomic.dt1.py and Borice.genomic.dt2.py).  This package should replace that one. The primary changes for version 3:

1.  The model estimates a separate value for t (the outcrossing rate) for each population specified in the popfile.  The v3 program outputs a posterior density for each population.
2.  The v3 program allows the user to run the program without data to investigate the prior distributions implied by parameter settings specified in Control.txt
3.  The v3 program takes a replication ID as a command line argument.  This facilitates parallel use of the program (different chains run simultaneously with outputs to files specific to the replication ID).
4.  The v3 program takes a standardized input for the genotype data.  The data for each SNP is a listed sequentially in the [name].genotypes.txt file.  Here, we include "Macrov2.genotypes.txt" as an example.  The format is exactly repeated for each SNP in the dataset:   

SNP	0	A	G

1	par	1.0	1.0	1.0

1	off	0.0	1.0	0.0

1	off	1.0	1.0	1.0

1	off	1.0	0.0825776470588	0.0

.
.
.

The first line of a "snp block" is the SNP ID with alternative bases listed.  Next, each familiy is listed in sequence (family ID is the first column).  With each family, the maternal plant ("par" in second column) is reported first, followed by each offspring in that family ("off" in second column).  The genotype likelihoods are reported in columns 3-5 for RR, RA, and AA, respectively.  Individuals with missing data have likelihood values of 1.0 1.0 1.0.  Maternal individuals are always included even if they are not genotypes.  In this case, they always have 1.0 1.0 1.0.

5.  The v3 program provides an additional output for inter-population divergence (Dxy for each pairwise contrast of populations).


Borice.Genomic.v3.py written in Python 2.7 (http://www.python.org/) and tested using the Linux operating system.  The program is dependent on SciPy.
If libraries are accessible, the programs should run using Python 2.7 in Windows or Mac operating systems. The command line usage is:

python Borice.Genomic.v3.py

The user needs to include two addtional files in the working directory: [name].genotypes.txt and the popfile.  The user needs to indentify these files by editing "Control.v3.txt".  The user can also specify design parameters in this file.  The particular example of Control.v3.txt contained here is set to run Borice on the genotypes in "Macrov2.genotypes.txt" using the popfile "macro.pops.v2.txt".

----------------

The AYMMr.tar.gz is bundle containing the AYMMr.py SNP filtering program in addition to two sample input files.  The command line usage for these test inputs is:

python AYMMr.py GL.example.txt key.example.txt

That run should produce an output file "AYMMr.scores.txt", a text file with the evaluation of each SNP in the input GL file.  The columns of the output: 

scaffold

SNP position

Number mother-offspring pairs 

P sub 1

P sub 2

Delta

Number of (not mother)-offspring pairs 

P sub 1 Not Mom

P sub 2 Not Mom

Delta Not Mom

-------------------------------

The second processing script, locus.police.py, is contained in the bundle "locus_police.tar.gz" along with sample input and output files.  Usage of the program for the sample inputs:

python locus.police.py hardcall

This should produce two output files: hardcall.impossibles.byMomOff.txt and hardcall.impossibles.bySNP.txt.  The SNP text file is a summary of impossible mother-offspring combinations per SNP.  Users my choose to simply delete all SNPs that exhibit incompatibilities prior to running BORICE or adopt an alternative calling method for Genotype Likelihoods.  The second output enumerates the number of incompatibilities per mother-offspring pairing across all SNPs.  Particular pairings that yield many incompatibilities may incorrectly assigned as mother and offspring.


--------------------------------------

Population.test.tar.gz contains the program to test for differences in outcrossing rate among populations (Population.t.py) as well as a sample input and output file.  Usage to run on the sample data:

python Population.t.py t.by.pop.txt

The input file is made from the out/self posterior probabilities from BORICE.



