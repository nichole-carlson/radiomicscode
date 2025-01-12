This repo contains the code for generating the analysis of the GRADS cohort.

GRADS imaging and clinical data need to be directly requested from the coordinating center at U Pittsburgh by emailing stevewis@pitt.edu, the head of the GRADS coordinating center.
0a_Clean_GRADSclinical.R - contains code to read in and select the demographic, PFT, PRO data on each individual.
0b_Clean_GRADSVAS.R - contains code to read in and select the visual assessment data.
0c_BuildComprehensiveGRADSData.R - merges the GRADS files prior to the radiomics.
0d_radiomicscreation.R - creates a radiomics dataset from pre-processed imaging data. Processing includes having a lung mask.
0e_MergeData.R - builds full analysis data.
0f_Filter_and_Cluster.R - contains a script to get cluster labels for the radiomics panel based on both using a decorrelation filter and then using RSKC means. To recreate the dataset in the shared drive, run this script.
1_Tables_1_and_2.R - contains code for tables with demographic information
2_Fig2and3_outcome_forest_plots.R - has code for creating forest plots of the regression coefficients by outcome for each radiomic grouping
3_Figure_1.R - makes a heatmap for the radiomic panel4_ScatterFigure.R - investigates the separation for the radiomics clustering
5_radiomics_var_Importance_analysis.R - code to identify the importance of the variables creating the clusters.
6_Sensitivity_Analysis1.R 
6_Sensitivity_Analysis2.R
7_Lung_Image_Figure_helperfunctions.R - functions to help with the figure of the different example image files.
7_Lung_Image_Figure.R - code that makes the panel of 9 lungs showing differences in 3 of our "important" variables
8_etable3_Oberstein.R - Supplemental table code.
9_summary_df_ibsi.r - Code to create the supplemental table on imaging parameters.
f_variable_names.R

Any custom functions used within these scripts can be found in the functions folder.