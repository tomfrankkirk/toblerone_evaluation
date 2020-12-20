Analysis scripts for the paper 'Toblerone: surface-based partial volume estimation'

There are three directories, named after the experiments they contain respectively, and one top-level script, run_evaluation.py which is the main analysis script. This script will run other scripts within the directories and produce 3 .mat files containing the results from each experiment. Finally, if you run the respective *_analysis.py scripts in each directory with the VS Code for data science plugin, it will produce and save the figures presented within the paper. 

For reasons of size, the intermediate files are not included here. Instead, the minimal set of files required to produce the intermediates are given, with the exception of HCP data which may be obtained directly from the HCP itself. In the top level script (run_evaluation), a variable  called ROOT is defined. This must point to a directory containing three subdirs, entitled sim_surfaces, brainweb and HCP_retest, within which the raw files for the experiments are located. 


