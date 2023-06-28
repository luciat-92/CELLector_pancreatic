# CELLector_pancreatic
Application of CELLector to pancreatic cancer (somatic mutations) 
and projection onto cell lines from CMAP and organoids lines from Corbo's lab

0. CELLector_preliminary.R: exploratory, can be skipped
1. CELLector_search_and_projection.R: search space from combined primary tumours and projection onto cell lines and organoid lines
2. convert_gdsctool_format.R: convert cell_ector output in proper format to perform pharmacogenomic analysis
3. gdsctools_drug_sens.py: drug sensitivity from gdsctools python library (ANOVA)
4. plot_drug_sensitivity_analysis.R: drug sensitivity boxplots and volcano plot
5. clinical_analysis_primary_tumors.R: patient survival analysis with Kaplan-Meier plots
