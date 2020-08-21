# MotionStructuralConnectivity

This repository contains all the code to reproduce the results and figures in Oldham et al., 2020. The efficacy of different preprocessing steps in reducing motion-related confounds in diffusion MRI connectomics. The paper can be found [here](https://www.sciencedirect.com/science/article/pii/S1053811920307382).

MAIN_ANALYSIS.m will perform all of the analyses done in the paper. 

PLOT_FIGURES.m will reproduce all the plots in the figure.

Data can be found [here](https://doi.org/10.26180/5e7313d012cee). Unzip the directories into the main directory so the scripts can locate them.

The exact code we used to generate all the tractograms/network matrices is contained in the QCSC_preprocessing.sh script. Note this script is mostly setup to work with how our data was arranged, so if you want to use it yourself it will require some tweaking.

Any issues please email stuart.oldham@monash.edu
