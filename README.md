# MotionStructuralConnectivity

This repository contains all the code to reproduce the results and figures in Oldham et al., 2020. The efficacy of different preprocessing steps in reducing motion-related confounds in diffusion MRI connectomics. A preprint of the paper can be found [here](https://www.biorxiv.org/content/10.1101/2020.03.25.008979v1).

MAIN_ANALYSIS.m will perform all of the analyses done in the paper. 

PLOT_FIGURES.m will reproduce all the plots in the figure.

Data can be found [here](https://figshare.com/s/3310385f29a156c93ca3). Unzip the directories into the main directory so the scripts can locate them.

The exact code we used to generate all the tractograms/network matrices is contained in the QCSC_preprocessing.sh script. Note this script is mostly setup to work with how our data was arranged, so if you want to use it yourself it will require some tweaking.

Any issues please email stuart.oldham@monash.edu
