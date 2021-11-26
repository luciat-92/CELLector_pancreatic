from gdsctools import *
import os

def getfile(filename, where='../data/'):
    return os.sep.join([where, filename])

def anova_report_save(ic50_file, gf_file, drug_file, out_file):
    
    gdsc = ANOVA(getfile(ic50_file), getfile(gf_file),
                 getfile(drug_file))

    results = gdsc.anova_all()
    results.to_csv(out_file)

    return results

# CS signatures
res_CS = anova_report_save(ic50_file = 'GDSC2_IC50.csv', 
                           gf_file = 'GDSC2_GF_CS_cmap_cl.csv', 
                           drug_file = 'GDSC2_drug_annotation.csv', 
                           out_file = '../output/gdsctools_CS_cmap_cl_output.csv')

# CFE signatures
res_CFE = anova_report_save(ic50_file = 'GDSC2_IC50.csv',
                            gf_file = 'GDSC2_GF_CFE_cmap_cl.csv', 
                            drug_file = 'GDSC2_drug_annotation.csv', 
                            out_file = '../output/gdsctools_CFE_cmap_cl_output.csv')             

