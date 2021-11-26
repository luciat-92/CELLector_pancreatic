# CELLector_pancreatic
Application of CELLector to pancreatic cancer (somatic mutations) 
and projection onto cell lines from CMAP and organoids lines from Corbo's lab

- CELLector_preliminary.R: exploratory, can be skipped
- CELLector_search_and_projection.R: search space from combined primary tumors and projection onto cell lines and organoid lines
- convert_gdsctool_format.R: convert cell_ector ouptu in proper format to perform pharmacogenomic analysis
- gdsctools_drug_sens.py: drug sensitivity from gdsctools python library (ANOVA)
- plot_drug_sensitivity_analysis.R: drug sensitivity boxplots and volcano plot
- clinical_analysis_primary_tumors.R: patient survival analysis with Kaplan-Meier plots