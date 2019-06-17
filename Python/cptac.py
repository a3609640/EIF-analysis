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

gene = 'EIF4E'
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

###############################################################################
# Use Case 2:  Comparing Clinical Data in endometrial cancer
### load the dataframe for clinical results by calling the en.get_clinical() method 
clinical_data = en.get_clinical()
print(clinical_data.columns)

# Use Case 3: Associating Clinical Variables with Proteomics
## Choose Clinical Attribute and Merge Dataframes
clinical_attribute = "tumor_Stage-Pathological"


## Merge attribute with roteomics dataframe
clinical_and_proteomics = en.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "proteomics", 
        metadata_cols    = clinical_attribute)
clinical_and_proteomics.head()

## Show possible variations of Histologic_type
clinical_and_proteomics[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_proteomics.filter(like='EIF4').columns)
#Index(['ANKHD1-EIF4EBP3_proteomics', 'EIF4A1_proteomics', 'EIF4A2_proteomics',
#       'EIF4A3_proteomics', 'EIF4B_proteomics', 'EIF4E_proteomics',
#       'EIF4E2_proteomics', 'EIF4E3_proteomics', 'EIF4EBP1_proteomics',
#       'EIF4EBP2_proteomics', 'EIF4EBP3_proteomics', 'EIF4ENIF1_proteomics',
#       'EIF4G1_proteomics', 'EIF4G2_proteomics', 'EIF4G3_proteomics',
#       'EIF4H_proteomics'],
#      dtype='object')

# plot graph on EIF4 proteomics
graphingSite = 'EIF4E_proteomics'
graphingSite = 'EIF4A1_proteomics'
graphingSite = 'EIF4G1_proteomics'
graphingSite = 'EIF4EBP1_proteomics'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_proteomics, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_proteomics, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('endometrial cancer')
#color = '.3' makes the dots black


## Merge attribute with phosphorylation dataframe
clinical_and_phosphorylation = en.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "phosphoproteomics", 
        metadata_cols    = clinical_attribute)
clinical_and_phosphorylation.head()

## Show possible variations of Histologic_type
clinical_and_phosphorylation[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_phosphorylation.filter(like='EIF4E').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4E-S24_phosphoproteomics'
graphingSite = 'EIF4EBP1-T36_phosphoproteomics'
graphingSite = 'EIF4EBP1-T70_phosphoproteomics'
graphingSite = 'EIF4EBP1-T46_phosphoproteomics'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_phosphorylation, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_phosphorylation, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('endometrial cancer')
#color = '.3' makes the dots black

## Merge attribute with phosphoproteomics_gene dataframe
clinical_and_phosphoproteomics_gene = en.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "phosphoproteomics_gene", 
        metadata_cols    = clinical_attribute)
clinical_and_phosphoproteomics_gene.head()

## Show possible variations of Histologic_type
clinical_and_phosphoproteomics_gene[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_phosphoproteomics_gene.filter(like='EIF4A').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4E_phosphoproteomics_gene'
graphingSite = 'EIF4E2_phosphoproteomics_gene'
graphingSite = 'EIF4G1_phosphoproteomics_gene'
graphingSite = 'EIF4EBP1_phosphoproteomics_gene'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_phosphoproteomics_gene, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_phosphoproteomics_gene, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('endometrial cancer')
#color = '.3' makes the dots black








###############################################################################
# Use Case 2:  Comparing Clinical Data in colon cancer
### load the dataframe for clinical results by calling the en.get_clinical() method 
clinical_data = col.get_clinical()
print(clinical_data.columns)

# Use Case 3: Associating Clinical Variables with Proteomics
## Choose Clinical Attribute and Merge Dataframes
clinical_attribute = "Stage"

## Merge attribute with roteomics dataframe
clinical_and_proteomics = col.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "proteomics", 
        metadata_cols    = clinical_attribute)
clinical_and_proteomics.head()

## Show possible variations of Histologic_type
clinical_and_proteomics[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_proteomics.filter(like='EIF4').columns)
#Index(['ANKHD1-EIF4EBP3_proteomics', 'EIF4A1_proteomics', 'EIF4A2_proteomics',
#       'EIF4A3_proteomics', 'EIF4B_proteomics', 'EIF4E_proteomics',
#       'EIF4E2_proteomics', 'EIF4E3_proteomics', 'EIF4EBP1_proteomics',
#       'EIF4EBP2_proteomics', 'EIF4EBP3_proteomics', 'EIF4ENIF1_proteomics',
#       'EIF4G1_proteomics', 'EIF4G2_proteomics', 'EIF4G3_proteomics',
#       'EIF4H_proteomics'],
#      dtype='object')

# plot graph on EIF4 proteomics
graphingSite = 'EIF4E_proteomics'
graphingSite = 'EIF4A1_proteomics'
graphingSite = 'EIF4G1_proteomics'
graphingSite = 'EIF4EBP1_proteomics'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_proteomics, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_proteomics, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('colon cancer')
#color = '.3' makes the dots black


## Merge attribute with phosphorylation dataframe
clinical_and_phosphorylation = col.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "phosphoproteomics", 
        metadata_cols    = clinical_attribute)
clinical_and_phosphorylation.head()

## Show possible variations of Histologic_type
clinical_and_phosphorylation[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_phosphorylation.filter(like='EIF4E').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4EBP1_T37__Q13541_phosphoproteomics'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_phosphorylation, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_phosphorylation, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('colon cancer')
#color = '.3' makes the dots black



clinical_and_transcriptomics = col.append_metadata_to_omics(
        metadata_df_name = "clinical", 
        omics_df_name    = "transcriptomics", 
        metadata_cols    = clinical_attribute)
clinical_and_transcriptomics.head()

## Show possible variations of Histologic_type
clinical_and_transcriptomics[clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(clinical_and_transcriptomics.filter(like='EIF4').columns)
#Index(['ANKHD1-EIF4EBP3_transcriptomics', 'EIF4A1_transcriptomics',
#       'EIF4A2_transcriptomics', 'EIF4A3_transcriptomics',
#       'EIF4B_transcriptomics', 'EIF4E_transcriptomics',
#       'EIF4E2_transcriptomics', 'EIF4E3_transcriptomics',
#       'EIF4EBP1_transcriptomics', 'EIF4EBP2_transcriptomics',
#       'EIF4EBP3_transcriptomics', 'EIF4ENIF1_transcriptomics',
#       'EIF4G1_transcriptomics', 'EIF4G2_transcriptomics',
#       'EIF4G3_transcriptomics', 'EIF4H_transcriptomics'],
#      dtype='object')

# plot graph on EIF4 transcriptomics
graphingSite = 'EIF4E_transcriptomics'
graphingSite = 'EIF4A1_transcriptomics'
graphingSite = 'EIF4G1_transcriptomics'
graphingSite = 'EIF4EBP1_transcriptomics'
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_transcriptomics, 
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_transcriptomics, 
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('colon cancer')
#color = '.3' makes the dots black















## Make dataframes with only endometrioid and only serous data in order to compare 
endom = clinical_and_phosphorylation.loc[
        clinical_and_phosphorylation[clinical_attribute] == "Endometrioid"]
serous = clinical_and_phosphorylation.loc[
        clinical_and_phosphorylation[clinical_attribute] == "Serous"]
## Here is where we set the NaN values to "Non_Tumor"
clinical_and_phosphorylation[[clinical_attribute]] = clinical_and_phosphorylation[[clinical_attribute]].fillna(value = "Non_Tumor")
## Make dataframes with non-tumor data as well
non = clinical_and_phosphorylation.loc[
        clinical_and_phosphorylation[clinical_attribute] == "Non_Tumor"]


num_sites_sufficient_data = 0
for num in range(1,len(serous.columns)):
    site = serous.columns[num]
    oneSite = serous[site]
    num_datapoints = oneSite.count() #Here we count the number of datapoints that we have at each site
    if num_datapoints >= 5:
        num_sites_sufficient_data += 1
print("Number of sites with sufficient data: ", num_sites_sufficient_data)
threshold = .05 / num_sites_sufficient_data
print("Multiple hypothesis adjusted threshold: ", threshold)


# find a site that is significantly different between non-tumor and serous
significantTests = []
significantSites = []
sigSiteCount = 0
p_value_index = 1
np.warnings.filterwarnings('ignore')
for num in range(1,len(clinical_and_phosphorylation.columns)):
    site = clinical_and_phosphorylation.columns[num]
    ttestRes = scipy.stats.ttest_ind(non[site], 
                                     serous[site])#This returns a tuple, with the t-statistic as 
    #the first value and the p-value as the second(index 1). We are interested in the p-value.
    if (ttestRes[p_value_index] <= threshold): #Check if there is a significant enough difference between data points
        sigSiteCount += 1
        significantTests.append(ttestRes)
        significantSites.append(site)
print("Number of sites with statistically significant difference between endometrioid and serous values: " + str(sigSiteCount))
print(significantSites)










# Find column whose name contains a EIF4E
print(clinical_and_phosphorylation.filter(like='EIF4E').columns)

# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4EBP1-T36_phosphoproteomics'
print(scipy.stats.ttest_ind(non[graphingSite], 
                            serous[graphingSite]))
sns.boxplot(x          = clinical_attribute, 
            y          = graphingSite, 
            data       = clinical_and_phosphorylation, 
            showfliers = False)
sns.stripplot(x        = clinical_attribute, 
              y        = graphingSite, 
              data     = clinical_and_phosphorylation, 
              color    = '.3')
#color = '.3' makes the dots black


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





