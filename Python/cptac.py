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

################################################################################
# Use Case 1: Correlation between Omics Data (Proteomics vs Transcriptom)
################################################################################
proteomics = en.get_proteomics()
samples = proteomics.index
proteins = proteomics.columns
print("Samples:", samples[0:20].tolist()) #the first twenty samples
print("Proteins:", proteins[0:20].tolist()) #the first twenty proteins
proteomics.head()

transcriptomics = en.get_transcriptomics()
transcriptomics.head()

## Step 1: Merging dataframes
SCD_cross = en.compare_omics(omics_df1_name = "proteomics",
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
         title  = 'Proteomics vs. Transcriptomics for ' + gene + ' in Endometrial Cancer')
plt.show()


################################################################################
# Use Case 2: Associating Clinical Variables with Proteomics in Endometrial Cancer
################################################################################
## load the dataframe for clinical results by calling the en.get_clinical() method
en_clinical_data = en.get_clinical()
print(en_clinical_data.columns)
## Choose clinical attribute and merge dataframes
en_clinical_attribute = "tumor_Stage-Pathological"


## Merge clinical attribute with proteomics dataframe
en_clinical_and_proteomics = en.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "proteomics",
        metadata_cols    = en_clinical_attribute)
en_clinical_and_proteomics[en_clinical_attribute] = en_clinical_and_proteomics[en_clinical_attribute].fillna("Normal")
en_clinical_and_proteomics.head()

## Show possible variations of Histologic_type
en_clinical_and_proteomics[en_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(en_clinical_and_proteomics.filter(like='JUN').columns)
#Index(['ANKHD1-EIF4EBP3_proteomics', 'EIF4A1_proteomics', 'EIF4A2_proteomics',
#       'EIF4A3_proteomics', 'EIF4B_proteomics', 'EIF4E_proteomics',
#       'EIF4E2_proteomics', 'EIF4E3_proteomics', 'EIF4EBP1_proteomics',
#       'EIF4EBP2_proteomics', 'EIF4EBP3_proteomics', 'EIF4ENIF1_proteomics',
#       'EIF4G1_proteomics', 'EIF4G2_proteomics', 'EIF4G3_proteomics',
#       'EIF4H_proteomics'],
#      dtype='object')

# plot graph on EIF4 proteomics
graphingSite = 'EIF4EBP1_proteomics'
graphingSite = 'EIF4A1_proteomics'
graphingSite = 'EIF4G1_proteomics'
graphingSite = 'JUN_proteomics'
sns.set_style("white")
sns.boxplot(x          = en_clinical_attribute,
            y          = graphingSite,
            data       = en_clinical_and_proteomics,
            showfliers = False,
           # order      = ['Tumor','Adjacent_normal','Myometrium_normal','Enriched_normal'],
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = en_clinical_attribute,
              y        = graphingSite,
              data     = en_clinical_and_proteomics,
              color    = '.3',
             # order    = ['Tumor','Adjacent_normal','Myometrium_normal','Enriched_normal'],
              order      = ["Normal","Stage I", "Stage II", "Stage III", "Stage IV"])
plt.xticks(rotation = 45)
plt.title('endometrial cancer')
#color = '.3' makes the dots black


## Merge clinical attribute with phosphorylation dataframe
en_clinical_and_phosphorylation = en.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "phosphoproteomics",
        metadata_cols    = en_clinical_attribute)
en_clinical_and_phosphorylation[en_clinical_attribute] = en_clinical_and_phosphorylation[en_clinical_attribute].fillna("Normal")
en_clinical_and_phosphorylation.head()

## Show possible variations of Histologic_type
en_clinical_and_phosphorylation[en_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(en_clinical_and_phosphorylation.filter(like='JUN-').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4E-S24_phosphoproteomics'
graphingSite = 'EIF4EBP1-T37_phosphoproteomics'
graphingSite = 'EIF4EBP1-T70_phosphoproteomics'
graphingSite = 'JUN-T93_phosphoproteomics'

sns.set_style("white")
sns.boxplot(x          = en_clinical_attribute,
            y          = graphingSite,
            data       = en_clinical_and_phosphorylation,
            showfliers = False,
            order      = ["Normal","Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = en_clinical_attribute,
              y        = graphingSite,
              data     = en_clinical_and_phosphorylation,
              color    = '.3',
              order    = ["Normal","Stage I", "Stage II", "Stage III", "Stage IV"])
plt.xticks(rotation = 45)
plt.title('endometrial cancer')
#color = '.3' makes the dots black


## Merge attribute with phosphoproteomics_gene dataframe
en_clinical_and_phosphoproteomics_gene = en.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "phosphoproteomics_gene",
        metadata_cols    = en_clinical_attribute)
en_clinical_and_phosphoproteomics_gene[en_clinical_attribute] = en_clinical_and_phosphoproteomics_gene[en_clinical_attribute].fillna("Normal")
en_clinical_and_phosphoproteomics_gene.head()

## Show possible variations of Histologic_type
en_clinical_and_phosphoproteomics_gene[en_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(en_clinical_and_phosphoproteomics_gene.filter(like='EIF4A').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4E_phosphoproteomics_gene'
graphingSite = 'JUN_phosphoproteomics_gene'
graphingSite = 'EIF4G1_phosphoproteomics_gene'
graphingSite = 'EIF4EBP1_phosphoproteomics_gene'

sns.set_style("white")
sns.boxplot(x          = en_clinical_attribute,
            y          = graphingSite,
            data       = en_clinical_and_phosphoproteomics_gene,
            showfliers = False,
            order      = ["Normal","Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = en_clinical_attribute,
              y        = graphingSite,
              data     = en_clinical_and_phosphoproteomics_gene,
              color    = '.3',
              order    = ["Normal","Stage I", "Stage II", "Stage III", "Stage IV"])
plt.xticks(rotation = 45)
plt.title('endometrial cancer')
#color = '.3' makes the dots black


################################################################################
# Use Case 3:  Comparing Clinical Data in colon cancer
################################################################################
### load the dataframe for clinical results by calling the en.get_clinical() method
col_clinical_data = col.get_clinical()
print(col_clinical_data.columns)

####################################################
## Associating Clinical Variables with proteomics ##
####################################################
## Choose Clinical Attribute and Merge Dataframes
col_clinical_attribute = "Stage"

## Merge clinical attribute with proteomics dataframe
col_clinical_and_proteomics = col.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "proteomics",
        metadata_cols    = col_clinical_attribute)


col_clinical_and_proteomics[col_clinical_attribute] = col_clinical_and_proteomics[col_clinical_attribute].fillna("Normal")

col_clinical_and_proteomics.head()

## Show possible variations of Histologic_type
col_clinical_and_proteomics[col_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(col_clinical_and_proteomics.filter(like='RAF').columns)
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
graphingSite = 'BRAF_proteomics'

sns.set_style("white")
sns.boxplot(x          = col_clinical_attribute,
            y          = graphingSite,
            data       = col_clinical_and_proteomics,
            showfliers = False,
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = col_clinical_attribute,
              y        = graphingSite,
              data     = col_clinical_and_proteomics,
              color    = '.3',
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"])
plt.xticks(rotation = 45)
plt.title('colon cancer')
#color = '.3' makes the dots black

###########################################################
## Associating Clinical Variables with Phosphoproteomics ##
###########################################################
## Merge attribute with phosphorylation dataframe
col_clinical_and_phosphorylation = col.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "phosphoproteomics",
        metadata_cols    = col_clinical_attribute)
col_clinical_and_phosphorylation[col_clinical_attribute] = col_clinical_and_phosphorylation[col_clinical_attribute].fillna("Normal")
col_clinical_and_phosphorylation.head()

## Show possible variations of Histologic_type
col_clinical_and_phosphorylation[col_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(col_clinical_and_phosphorylation.filter(like='PIK3C').columns)
# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4EBP1_S101__Q13541_phosphoproteomics'

sns.set_style("white")
sns.boxplot(x          = col_clinical_attribute,
            y          = graphingSite,
            data       = col_clinical_and_phosphorylation,
            showfliers = False,
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = col_clinical_attribute,
              y        = graphingSite,
              data     = col_clinical_and_phosphorylation,
              color    = '.3',
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"])
plt.xticks(rotation = 45)
plt.title('colon cancer')
#color = '.3' makes the dots black



col_clinical_and_transcriptomics = col.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "transcriptomics",
        metadata_cols    = clinical_attribute)
col_clinical_and_transcriptomics.head()

## Show possible variations of Histologic_type
col_clinical_and_transcriptomics[col_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(col_clinical_and_transcriptomics.filter(like='EIF4').columns)
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
sns.boxplot(x          = col_clinical_attribute,
            y          = graphingSite,
            data       = col_clinical_and_transcriptomics,
            showfliers = False,
            order      = ["Stage I", "Stage II", "Stage III", "Stage IV"])
sns.stripplot(x        = col_clinical_attribute,
              y        = graphingSite,
              data     = col_clinical_and_transcriptomics,
              color    = '.3',
              order    = ["Stage I", "Stage II", "Stage III", "Stage IV"])
plt.title('colon cancer')
#color = '.3' makes the dots black


##########################################################################
## Associating Clinical Variables with proteomics and Phosphoproteomics ##
##########################################################################
col_clinical_and_proteomics_phosphoproteomics = pd.concat(
        [col_clinical_and_proteomics, col_clinical_and_phosphorylation], 
        axis=1, sort=False) 














################################################################################
# Use Case 4:  Comparing Clinical Data in ovarian cancer
################################################################################
### load the dataframe for clinical results by calling the en.get_clinical() method
ov_clinical_data = ov.get_clinical()
print(ov_clinical_data.columns)

####################################################
## Associating Clinical Variables with Proteomics ##
####################################################
## Choose Clinical Attribute and Merge Dataframes
ov_clinical_attribute = "Tumor_Stage_Ovary_FIGO"

## Merge clinical attribute with proteomics dataframe
ov_clinical_and_proteomics = ov.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "proteomics",
        metadata_cols    = ov_clinical_attribute)
ov_clinical_and_proteomics[ov_clinical_attribute] = ov_clinical_and_proteomics[ov_clinical_attribute].fillna("Normal")
ov_clinical_and_proteomics.head()

## Show possible variations of Histologic_type
ov_clinical_and_proteomics[ov_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(ov_clinical_and_proteomics.filter(like='EIF4G1').columns)
#Index(['ANKHD1-EIF4EBP3_proteomics', 'EIF4A1_proteomics', 'EIF4A2_proteomics',
#       'EIF4A3_proteomics', 'EIF4B_proteomics', 'EIF4E_proteomics',
#       'EIF4E2_proteomics', 'EIF4E3_proteomics', 'EIF4EBP1_proteomics',
#       'EIF4EBP2_proteomics', 'EIF4EBP3_proteomics', 'EIF4ENIF1_proteomics',
#       'EIF4G1_proteomics', 'EIF4G2_proteomics', 'EIF4G3_proteomics',
#       'EIF4H_proteomics'],
#      dtype='object')
cols = []
count = 1
for column in ov_clinical_and_proteomics.columns:
    if column == 'EIF4G1_proteomics':
        cols.append('EIF4G1_proteomics'+ str(count))
        count+=1
        continue
    cols.append(column)
ov_clinical_and_proteomics.columns = cols

print(ov_clinical_and_proteomics.filter(like='JUN').columns)
cols = []
count = 1
for column in ov_clinical_and_proteomics.columns:
    if column == 'JUN_proteomics':
        cols.append('JUN_proteomics'+ str(count))
        count+=1
        continue
    cols.append(column)
ov_clinical_and_proteomics.columns = cols

# plot graph on EIF4 proteomics
graphingSite = 'EIF4E_proteomics'
graphingSite = 'EIF4A1_proteomics'
graphingSite = 'EIF4EBP1_proteomics'
graphingSite = 'JUN_proteomics1'

sns.set_style("white")
sns.boxplot(x          = ov_clinical_attribute,
            y          = graphingSite,
            data       = ov_clinical_and_proteomics,
            showfliers = False,
            order      = ["Normal", "IC", "IIIA", "IIIB", "IIIC","IV"])
sns.stripplot(x        = ov_clinical_attribute,
              y        = graphingSite,
              data     = ov_clinical_and_proteomics,
              color    = '.3',
            order      = ["Normal", "IC", "IIIA", "IIIB", "IIIC","IV"])
plt.xticks(rotation = 45)
plt.title('ovarian cancer')
#color = '.3' makes the dots black

###########################################################
## Associating Clinical Variables with Phosphoproteomics ##
###########################################################
## Merge attribute with phosphorylation dataframe
ov_clinical_and_phosphorylation = ov.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "phosphoproteomics",
        metadata_cols    = ov_clinical_attribute)
ov_clinical_and_phosphorylation[ov_clinical_attribute] = ov_clinical_and_phosphorylation[ov_clinical_attribute].fillna("Normal")
ov_clinical_and_phosphorylation.head()

## Show possible variations of Histologic_type
ov_clinical_and_phosphorylation[ov_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(ov_clinical_and_phosphorylation.filter(like='EIF4E-').columns)
cols = []
count = 1
for column in ov_clinical_and_phosphorylation.columns:
    if column == 'EIF4E-S24s_phosphoproteomics':
        cols.append('EIF4E-S24s_phosphoproteomics'+ str(count))
        count+=1
        continue
    cols.append(column)
ov_clinical_and_phosphorylation.columns = cols

print(ov_clinical_and_phosphorylation.filter(like='EIF4EBP1').columns)
cols = []
count = 1
for column in ov_clinical_and_phosphorylation.columns:
    if column == 'EIF4EBP1-T46t_phosphoproteomics':
        cols.append('EIF4EBP1-T46t_phosphoproteomics'+ str(count))
        count+=1
        continue
    cols.append(column)
ov_clinical_and_phosphorylation.columns = cols

cols = []
count = 1
for column in ov_clinical_and_phosphorylation.columns:
    if column == 'EIF4EBP1-T70t_phosphoproteomics':
        cols.append('EIF4EBP1-T70t_phosphoproteomics'+ str(count))
        count+=1
        continue
    cols.append(column)
ov_clinical_and_phosphorylation.columns = cols

# plot graph on EIF4 phosphorylation
graphingSite = 'EIF4EBP1-T37tT46t_phosphoproteomics'

sns.set_style("white")
sns.boxplot(x          = ov_clinical_attribute,
            y          = graphingSite,
            data       = ov_clinical_and_phosphorylation,
            showfliers = False,
            order      = ["Normal", "IC", "IIIA", "IIIB", "IIIC", "IV"])
sns.stripplot(x        = ov_clinical_attribute,
              y        = graphingSite,
              data     = ov_clinical_and_phosphorylation,
              color    = '.3',
            order      = ["Normal", "IC", "IIIA", "IIIB", "IIIC", "IV"])
plt.xticks(rotation = 45)
plt.title('ovarian cancer')
#color = '.3' makes the dots black


###########################################################
## Associating Clinical Variables with Transcriptomics ##
###########################################################
ov_clinical_and_transcriptomics = ov.append_metadata_to_omics(
        metadata_df_name = "clinical",
        omics_df_name    = "transcriptomics",
        metadata_cols    = ov_clinical_attribute)
ov_clinical_and_transcriptomics.head()

## Show possible variations of Histologic_type
ov_clinical_and_transcriptomics[ov_clinical_attribute].unique()

# Find column whose name contains a EIF4E
print(ov_clinical_and_transcriptomics.filter(like='EIF4G').columns)
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
sns.set_style("white")
sns.boxplot(x          = ov_clinical_attribute,
            y          = graphingSite,
            data       = ov_clinical_and_transcriptomics,
            showfliers = False,
            order      = ["Normal", "IC", "IIIA", "IIIB", "IIIC", "IV"])
sns.stripplot(x        = ov_clinical_attribute,
              y        = graphingSite,
              data     = ov_clinical_and_transcriptomics,
              color    = '.3',
              order    = ["Normal", "IC", "IIIA", "IIIB", "IIIC", "IV"])
plt.title('colon cancer')
#color = '.3' makes the dots black









