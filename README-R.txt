FaST-LMM-EWASher v1.0 README   2/14/2014

I. Content
----------
R/ 		R code for FaST-LMM-EWASher
demo/		Example dataset
README-R.txt	This document

II. Example
-----------
usage: FLE-script.r data_file phenotype_file [-covar covar_file] [-map map_file]

To run FaST-LMM-EWASher on the example dataset, run this command from the command line in directory .\demo:
rscript ..\R\FLE-script.r input_data.txt input_phenotype.txt -covar input_covariates.txt

To optionally specify the map file (see "Inputs" documentation below), use:
rscript ..\R\FLE-script.r input_data.txt input_phenotype.txt -covar input_covariates.txt -map input_testmarkers.map

III. Inputs
-----------
data_file
A tab-delimited table of methylation values. In 27k and 450k arrays, this correspond to the beta values of probes. Each column 
corresponds to a sample and each row corresponds to a probe. The first row is a header: "ID", sampleID1, sampleID2, ... Each 
subsequent row has the form: probeID, value1, value2, ...

phenotype_file
Each row has the form: sampleID	sampleStatus. Sample status is 1 if the sample is a case and is 0 otherwise. This file is 
tab-delimited. 

covar_file [optional]
The user may specify a set of covariates to include in the model (e.g. age, gender). Each row corresponds to a sample. The order
of the samples must be the same as in the data_file header. The first column is always set to '1', the second column is the 
sample ID, and each subsequent column correspond to a covariate. Do not add column header. This file is tab-delimited.

map_file [optional]
This program requires the chromosome and genetic distance of every marker. By default, it assumes a 450K/27K Illumina methylation
chip, in which case this parameter can be left blank. If the default is not applicable, please provide a map file, which contains
the [chrm#,name,genetic_distance] of each marker, oner per row, with no header, and optionally, extra columns to be ignored. 
An example is provided in the demo directory (testmarkers.map).

IV. Output
-----------
FaST-LMM-EWASher creates two output folders, results/ and tmp/.

The folder results/ contains the following files:
out_ewasher.txt. The EWAS P values and statistics for each loci as computed by FaST-LMM-EWASher.  FaST-LMM-EWASher filters out 
	loci that are constitutively on or off, and the association statistics are computed for the remaining loci.
out_linreg.txt. The P values and statistics from a linear regression with the covariates specified by user. 
qq_ewasher.png. QQ-plot of the FaST-LMM-EWASher P values.
qq_linreg.png. QQ-plot of the linear regression P values.
similarity.txt. The similarity matrix computed by FaST-LMM-EWASher.
summary.txt. A summary file that states how many loci where analyzed and how many principle components were used by 
	FaST-LMM-EWASher.

The folder tmp/ contains intermediate files produced by FaST-LMM-EWASher.

V. Debugging
-------------
The file tmp/log.txt contains more detailed run time reports. 

FaST-LMM-EWASher uses Intel's MKL library to perform math operations. The program gives an error message if the MKL on the 
computer is outdated. If this occurs, set the environment variable FastLmmUseAnyMklLib=1. For example, in Unix, add 
"export FastLmmUseAnyMklLib=1" in the .bashrc or .cshrc file.     

VI. Miscellaneous Notes
------------------------
By default, this software uses REML to estimate the linear mixed model parameters, in conjunction with an F-test to
obtain p-values.

There are a few parameters currently hard-coded, at the top of fastlmm-ewasher.py and can be changed there:
#constitutively unmethylated/methylated filter values
minFilterVal = 0.2; maxFilterVal = 0.8 

#values for fastlmmc autoselect to try for number of loci in GSM
#may be advantageous to try a finer, or more enlarged grid, but this default setting seems reasonable.
asVals = append(append(c(1), seq(100, 900, 100)), seq(1000, 10000, 1000))

# maximum number of pcs to allow
maxPC = 10

#threshold on lambda for when to stop adding PCs
lamdaThresh = 1.0
