README

This folder contains script files for the male-female infertility genetic interaction analysis.
The software used are PLINK and R.

Script descriptions:
MFchr.sh is the master shell script
MFinter.R is the main analysis script
MFmplot.R merges and plots results
All R scripts beginning with "f_" are functions

Input files required are binary genetic files (.bim, .bed, and .fam).
Also required is a key file linking parental couples for each pregnancy and specifying the phenotype (0, 1), and including a 'couple ID' to identify couples with multiple pregnancies.
A covariate file is optional.
Wherever an input file is expected in the code it is denoted by FILE.filetype, e.g. GENETICS.bim.
Paths are denoted by /PATH_DESCRIPTION/, e.g. /PATH_TO_SCRIPT_FILES/.

All code is written by Siri Nærland Skodvin (2022).
