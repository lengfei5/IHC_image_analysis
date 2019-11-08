## this file is to explain how to use the three scripts to segment images and to save results 
General:
There are two scripts: 
(1) "Image_processing.m": low-level function, meaning that you do not need to open it and leave them alone in the usual cases.
(2) "processing.pl":  the function you are going to use and to make configuration before using it.

Usage:/Volumes/DANIEL/IHC-OBESE/Restricted Feeding Exp Ob/Image_01.btf
You have to open your terminal and go to the folder where these scripts are. 
If you want to segment 5% of an image named "Image_01.tif" in the path of 
"/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images2run/", then type the following in the terminal:
perl -w processing.pl "/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images2run/Image_01.tif" 0.05

If you want to segment 15% of all images in the path of 
"/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images2run/", then type the following in the terminal:
perl -w processing.pl "/Users/jiwang/Proteomics_anaylysis/Nuclear_proteins/Liver_cell_nuclei/15_may_2015_IHC/images2run/*" 0.15

Configurations:
* make sure perl and Matlab installed and working
* put all images to analyze into one folder; create a new folder for the results
* modify the path to save the results in "processing.pl" 
	$dirout = "/USERS/JIWANG/PROTEOMICS_ANAYLYSIS/NUCLEAR_PROTEINS/LIVER_CELL_NUCLEI/15_MAY_2015_IHC/RESULTS/SAVED/";
* modify the path to Matlab in the "processing.pl"
	$command = "/Applications/MATLAB_R2014b.app/bin/matlab ..." 