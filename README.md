## **Scripts to prepare input files and pool compounds into mixtures.  Scripts have detailed information using the --help option.**
#### This workflow is broken down into 4 steps. You may not need all the steps.
##### Step1 
This script will take a file with plate and structure information and add an exact mass and molecular formula column.  This is calculated using RDKIT.

##### Step2 

This script is the main pooling script where it will take the output from Step1 which is a file containing plate information with mass and formula.  The script first sorts by exact mass and then sequentially places compounds into the plate.  This aims to give two results.  Pools with minimal mass clashes in each well, and a normally distributed mass values in each well.  The script also outputs diagnostic files to let the user know how many mass clashes are in each well as defined by the input mass resolution.  The user can adjust the number of compounds per well or number of wells to balance samples and mass clashes.  Performance depends on the distribution of mass values of the compounds as well as the resolution of the users mass spec analysis method. In this publication the resolution was >0.05 amu.

##### Step3 

This script converts the pooling output into a set of transfer instructions for Beckman Coulter and we can confirm usage on both 650 and 655 machines.

##### Step4 

This script converts the pooling output int a set of files which allows matching of hits to the compound structures.  The script writes a file for each or run which is analyzed. This format was designed to work with the Agilent Mass Hunter PCDL and Quantitative Analysis modules.  

