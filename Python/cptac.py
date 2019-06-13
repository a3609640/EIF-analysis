import cptac
import math
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
import seaborn as sns
import statsmodels.stats.multitest


# To view available datasets, enter 'cptac.list_data()'.
cptac.list_data()
en = cptac.Endometrial()
col = cptac.Colon()
ov = cptac.Ovarian()

en.list_data()
col.list_data()
ov.list_data()

proteomics = en.get_proteomics()
samples = proteomics.index
proteins = proteomics.columns
print("Samples:", samples[0:20].tolist()) #the first twenty samples
print("Proteins:", proteins[0:20].tolist()) #the first twenty proteins
proteomics.head()

transcriptomics = en.get_transcriptomics()
transcriptomics.head()



# Use Case 1: Comparing Omics Data
en = cptac.Endometrial()
## Step 1: Merging dataframes
SCD_cross = ov.compare_omics(omics_df1_name = "proteomics", 
                             omics_df2_name = "transcriptomics", 
                             genes1         = "SCD", 
                             genes2         = "SCD")
SCD_cross.head()
## Step 2: Plot data
sns.set(style = "darkgrid")
plot = sns.regplot(x    = SCD_cross.columns[0], 
                   y    = SCD_cross.columns[1], 
                   data = SCD_cross)
plot.set(xlabel = 'Proteomics', 
         ylabel = 'Transcriptomics', 
         title  = 'Proteomics vs. Transcriptomics for the SCD gene')
plt.show()

gene = 'EIF4G1'
gene_cross = col.compare_omics(omics_df1_name = "proteomics", 
                               omics_df2_name = "transcriptomics", 
                               genes1         = gene, 
                               genes2         = gene)
plot = sns.regplot(x     = gene_cross.columns[0], 
                   y     = gene_cross.columns[1], 
                   data  = gene_cross, 
                   color = "green")
plot.set(xlabel = 'Proteomics', 
         ylabel = 'Transcriptomics', 
         title  = 'Proteomics vs. Transcriptomics for the ' + gene + ' gene')
plt.show()


gene_cross = en.compare_omics(omics_df1_name = "proteomics", 
                              omics_df2_name = "transcriptomics", 
                              genes1         = gene, 
                              genes2         = gene)
plot = sns.regplot(x     = gene_cross.columns[0], 
                   y     = gene_cross.columns[1], 
                   data  = gene_cross, 
                   color = "green")
plot.set(xlabel = 'Proteomics', 
         ylabel = 'Transcriptomics', 
         title  = 'Proteomics vs. Transcriptomics for the ' + gene + ' gene')
plt.show()


gene_cross = ov.compare_omics(omics_df1_name = "proteomics", 
                              omics_df2_name = "transcriptomics", 
                              genes1         = gene, 
                              genes2         = gene)
plot = sns.regplot(x     = gene_cross.columns[0], 
                   y     = gene_cross.columns[1], 
                   data  = gene_cross, 
                   color = "green")
plot.set(xlabel = 'Proteomics', 
         ylabel = 'Transcriptomics', 
         title  = 'Proteomics vs. Transcriptomics for the ' + gene + ' gene')
plt.show()

# Use Case 3: Associating Clinical Variables with Acetylation













en = cptac.Endometrial()
somatic_mutations = en.get_mutations()

gene = somatic_mutations["Gene"].value_counts().index[0]
print(gene)

omics_mutations = en.append_mutations_to_omics(omics_df_name  = "proteomics", 
                                               mutation_genes = gene, 
                                               omics_genes    = gene)
omics_mutations.head()

a4_dims = (11.7, 8.27) #dimensions for bigger plot
fig, ax = plt.subplots(figsize = a4_dims) #bigger plot displays Somatic Gene Mutation category without overlapping labels
somatic_boxplot = sns.boxplot(data  = omics_mutations, 
                              x     = gene + "_Mutation",
                              y     = gene + "_proteomics", 
                              ax    = ax, 
                              order = 
                                      ["Wildtype_Normal",
                                       "Wildtype_Tumor",
                                       "Missense_Mutation",
                                       "Nonsense_Mutation",
                                       "Frame_Shift_Ins",
                                       "Frame_Shift_Del",
                                       "Splice_Site"]
                              ) #order parameter is used to reorder the mutation categories 
somatic_boxplot.set_title("PTEN gene mutation protein abundance")
ax = sns.boxplot(x = gene + "_Mutation", 
                 y = gene + "_proteomics", 
                 data = omics_mutations)


somatic_boxplot = sns.stripplot(data   = omics_mutations, 
                                x      = gene + "_Mutation",
                                y      = gene + "_proteomics",
                                jitter = True, 
                                color  = ".3", 
                                order  = 
                                        ["Wildtype_Normal",
                                         "Wildtype_Tumor",
                                         "Missense_Mutation",
                                         "Nonsense_Mutation",
                                         "Frame_Shift_Ins",
                                         "Frame_Shift_Del",
                                         "Splice_Site"]
                                )
somatic_boxplot.set(xlabel = "Somatic Gene Mutation",
                    ylabel = "Proteomics")
plt.show()


# Use Case 6: Comparing Derived Molecular Data with Proteomics
import pandas as pd
import seaborn as sns
import matplotlib.pyplot as plt
import cptac
en = cptac.Endometrial()





