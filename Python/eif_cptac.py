import cptac
import math
import matplotlib
import matplotlib.pyplot as plt
import pandas as pd
import numpy as np
import scipy.stats
import seaborn as sns
import statsmodels.stats.multitest
from statannot import add_stat_annotation
from scipy.stats import pearsonr


def downloadCptac():
    # To view available datasets, enter 'cptac.list_data()'.
    cptac.list_datasets()
    cptac.download(dataset = "endometrial")
    cptac.download(dataset = 'colon')
    cptac.download(dataset = 'ovarian')
    cptac.download(dataset = 'RenalCcrcc')
    #cptac.download(dataset ='luad')
    #cptac.download(dataset ='brca')
downloadCptac()

endometrialData = cptac.Endometrial()
colorectalData = cptac.Colon()
ovarianData = cptac.Ovarian()
renalData = cptac.RenalCcrcc()
lungData = cptac.Luad()
breastData = cptac.Brca()

def listDataForEachCancer():
    print("endometrial")
    endometrialData.list_data()
    print("\n\ncolorectal")
    colorectalData.list_data()
    print("\n\novarian")
    ovarianData.list_data()
    print("\n\nrenal")
    renalData.list_data()

listDataForEachCancer()

#################################################################
# Correlation: Proteomics vs Transcriptom in Endometrial Cancer #
#################################################################
def correlationPlot(dataSet, label, omics1, omics2, gene1, gene2):
    gene_cross = dataSet.join_omics_to_omics(df1_name = omics1,
                                             df2_name = omics2,
                                             genes1   = gene1,
                                             genes2   = gene2)
    print(gene_cross.head())
    gene_cross = gene_cross.dropna()
    corr = pearsonr(gene_cross.iloc[:, 0],
                    gene_cross.iloc[:, 1])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style      ="white",
            font_scale = 1.5)
    plot = sns.regplot(x    = gene_cross.columns[0],
                       y    = gene_cross.columns[1],
                       data = gene_cross)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize = 12)
    plot.set(xlabel = gene1 + ' ' + omics1,
             ylabel = gene2 + ' ' + omics2,
             title = '{} vs {} ({})'.format(gene1, gene2, label))
    plt.show()


correlationPlot(dataSet = endometrialData,
                label  = "Endometrial Cancer",
                omics1 = "proteomics",
                omics2 = "transcriptomics",
                gene1  = "EIF4A1",
                gene2  = "VEGFA")

correlationPlot(dataSet = endometrialData,
                label  = "Endometrial Cancer",
                omics1 = "phosphoproteomics_gene",
                omics2 = "transcriptomics",
                gene1  = "EIF4G1",
                gene2  = "EIF4G1")

correlationPlot(dataSet = endometrialData,
                label  = "Endometrial Cancer",
                omics1 = "phosphoproteomics_gene",
                omics2 = "proteomics",
                gene1  = "EIF4A1",
                gene2  = "EIF4A1")

correlationPlot(dataSet = colorectalData,
                label  = "Colorectal Cancer",
                omics1 = "proteomics",
                omics2 = "proteomics",
                gene1  = "EIF4G1",
                gene2  = "EIF4A1")

try:
    correlationPlot(dataSet = ovarianData,
                    label  = "Ovarian Cancer",
                    omics1 = "proteomics",
                    omics2 = "transcriptomics",
                    gene1  = "EIF4G1",
                    gene2  = "EIF4G1")
except Exception as ex:
    print('Could not make correlation plot for Ovarian Cancer: ' + str(ex))


######################################################################
# Correlation: Phosphoproteomics vs Proteomics in Endometrial Cancer #
######################################################################
## correlation between 4EBP1-T37 and 4EBP1 protein in endometrial cancer
def siteSpecificCorrelationPlot(dataSet, gene1, gene2, site):
    gene_cross_en = dataSet.join_omics_to_omics(
            df1_name = "phosphoproteomics",
            df2_name = "proteomics",
            genes1   = gene1,
            genes2   = gene2)
    print(gene_cross_en.columns)
    gene_cross_en = gene_cross_en.dropna(subset = ['EIF4EBP1-T37_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['EIF4EBP1-T37_phosphoproteomics'],
                    gene_cross_en[site])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style      ="white",
            font_scale = 1.5)
    plot = sns.regplot(x    = gene_cross_en['EIF4EBP1-T37_phosphoproteomics'],
                       y    = gene_cross_en[site],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'EIF4EBP1-T37 phosphoproteomics',
             ylabel = site,
             title  = 'Phosphoproteomics vs Proteomics (Endometrial Cancer)')
    plt.show()

siteSpecificCorrelationPlot(dataSet = endometrialData,
                  gene1 = "EIF4EBP1",
                  gene2 = "EIF4EBP1",
                  site = "EIF4EBP1_proteomics")

'''
## correlation between 4EBP1-T37 and MYC protein in endometrial cancer
def encorrelationplot(gene1, gene2):
    gene_cross_en = en.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "proteomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_en.head())
    gene_cross_en = gene_cross_en.dropna(subset = ['EIF4EBP1-T37_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['EIF4EBP1-T37_phosphoproteomics'],
                    gene_cross_en['EIF4A1_proteomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style      ="white",
            font_scale = 1.5)
    plot = sns.regplot(x    = gene_cross_en['EIF4EBP1-T37_phosphoproteomics'],
                       y    = gene_cross_en['EIF4A1_proteomics'],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'EIF4EBP1-T37 phosphoproteomics',
             ylabel = 'EIF4A1 proteomics',
             title  = 'Phosphoproteomics vs Proteomics (Endometrial Cancer)')
    plt.show()
'''
siteSpecificCorrelationPlot(dataSet = endometrialData,
                  gene1 = "EIF4EBP1",
                  gene2 = "EIF4A1",
                  site = "EIF4A1_proteomics")


#########################################################################
# Association: Clinical Variables with Proteomics in Endometrial Cancer #
#########################################################################

## load the dataframe for clinical results by calling the en.get_clinical() method
en_clinical_data = en.get_clinical()
print(en_clinical_data.columns)

## Clinical Variables with Proteomics
def en_pro_cli_plot(gene):
    en_clinical_and_proteomics = en.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "proteomics",
            metadata_cols    = "tumor_Stage-Pathological",
            omics_genes      = gene)
    en_clinical_and_proteomics["tumor_Stage-Pathological"] = en_clinical_and_proteomics["tumor_Stage-Pathological"].fillna("Normal")
    en_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    en_clinical_and_proteomics["tumor_Stage-Pathological"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
    ax = sns.boxplot(x          = "tumor_Stage-Pathological",
                     y          = gene + '_proteomics',
                     data       = en_clinical_and_proteomics,
                     showfliers = False,
                     order      = order)
    sns.stripplot(x        = "tumor_Stage-Pathological",
                  y        = gene + '_proteomics',
                  data     = en_clinical_and_proteomics,
                  color    = '.3',
                  order     = order)
    add_stat_annotation(ax,
                        data        = en_clinical_and_proteomics,
                        x           = "tumor_Stage-Pathological",
                        y           = gene + '_proteomics',
                        order       = order,
                        boxPairList = [("Normal", "Stage I"),
                                       ("Normal", "Stage II"),
                                       ("Normal", "Stage III"),
                                       ("Normal", "Stage IV")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('endometrial cancer')

en_pro_cli_plot(gene = "EIF4E")


## Clinical Variables with phosphoproteomics
def en_phos_cli_plot(gene):
    en_clinical_and_proteomics = en.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics",
            metadata_cols    = "tumor_Stage-Pathological",
            omics_genes      = gene)
    en_clinical_and_proteomics["tumor_Stage-Pathological"] = en_clinical_and_proteomics["tumor_Stage-Pathological"].fillna("Normal")
    en_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    en_clinical_and_proteomics["tumor_Stage-Pathological"].unique()
    PhosphoSite = list(en_clinical_and_proteomics.filter(like= gene).columns.values.tolist())
    for i in PhosphoSite:
        print(i)
        try:
            en_clinical_and_proteomics = en.join_metadata_to_omics(
                    metadata_df_name = "clinical",
                    omics_df_name    = "phosphoproteomics",
                    metadata_cols    = "tumor_Stage-Pathological",
                    omics_genes      = gene)
            en_clinical_and_proteomics["tumor_Stage-Pathological"] = en_clinical_and_proteomics["tumor_Stage-Pathological"].fillna("Normal")
            en_clinical_and_proteomics = en_clinical_and_proteomics.dropna(subset = [i])
            plt.figure()
            sns.set(style      ="white",
                    font_scale = 1.5)
            order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
            ax = sns.boxplot(x          = "tumor_Stage-Pathological",
                             y          = i,
                             data       = en_clinical_and_proteomics,
                             showfliers = False,
                             order      = order)
            sns.stripplot(x        = "tumor_Stage-Pathological",
                          y        = i,
                          data     = en_clinical_and_proteomics,
                          color    = '.3',
                          order     = order)
            add_stat_annotation(ax,
                                data        = en_clinical_and_proteomics,
                                x           = "tumor_Stage-Pathological",
                                y           = i,
                                order       = order,
                                boxPairList = [("Normal", "Stage I"),
                                               ("Normal", "Stage II"),
                                               ("Normal", "Stage III"),
                                               ("Normal", "Stage IV")],
                                test        = 't-test_ind',
                                textFormat  = 'star',
                                loc         = 'inside',
                                verbose     = 2)
            plt.title('endometrial cancer')
        except: ValueError
        pass

en_phos_cli_plot(gene = "PIK3CA")

## Clinical Variables with phosphoproteomics total
def en_phos_cli_plot(gene):
    en_clinical_and_proteomics = en.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics_gene",
            metadata_cols    = "tumor_Stage-Pathological",
            omics_genes      = gene)
    en_clinical_and_proteomics["tumor_Stage-Pathological"] = en_clinical_and_proteomics["tumor_Stage-Pathological"].fillna("Normal")
    en_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    en_clinical_and_proteomics["tumor_Stage-Pathological"].unique()
    PhosphoSite = list(en_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in PhosphoSite:
        print(i)
        en_clinical_and_proteomics = en_clinical_and_proteomics.dropna(subset = [i])
        plt.figure()
        sns.set(style      ="white",
                font_scale = 1.5)
        order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
        ax = sns.boxplot(x          = "tumor_Stage-Pathological",
                         y          = i,
                         data       = en_clinical_and_proteomics,
                         showfliers = False,
                         order      = order)
        sns.stripplot(x        = "tumor_Stage-Pathological",
                      y        = i,
                      data     = en_clinical_and_proteomics,
                      color    = '.3',
                      order     = order)
        add_stat_annotation(ax,
                            data        = en_clinical_and_proteomics,
                            x           = "tumor_Stage-Pathological",
                            y           = i,
                            order       = order,
                            boxPairList = [("Normal", "Stage I"),
                                           ("Normal", "Stage II"),
                                           ("Normal", "Stage III"),
                                           ("Normal", "Stage IV")],
                            test        = 't-test_ind',
                            textFormat  = 'star',
                            loc         = 'inside',
                            verbose     = 2)
        plt.title('endometrial cancer')

en_phos_cli_plot(gene = "EIF4EBP1")



## ## Clinical Variables with transcriptomics
def en_trans_cli_plot(gene):
    en_clinical_and_proteomics = en.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "transcriptomics",
            metadata_cols    = "tumor_Stage-Pathological",
            omics_genes      = gene)
    en_clinical_and_proteomics["tumor_Stage-Pathological"] = en_clinical_and_proteomics["tumor_Stage-Pathological"].fillna("Normal")
    en_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    en_clinical_and_proteomics["tumor_Stage-Pathological"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
    ax = sns.boxplot(x          = "tumor_Stage-Pathological",
                     y          = gene + '_transcriptomics',
                     data       = en_clinical_and_proteomics,
                     showfliers = False,
                     order      = order)
    sns.stripplot(x        = "tumor_Stage-Pathological",
                  y        = gene + '_transcriptomics',
                  data     = en_clinical_and_proteomics,
                  color    = '.3',
                  order     = order)
    add_stat_annotation(ax,
                        data        = en_clinical_and_proteomics,
                        x           = "tumor_Stage-Pathological",
                        y           = gene + '_transcriptomics',
                        order       = order,
                        boxPairList = [("Normal", "Stage I"),
                                       ("Normal", "Stage II"),
                                       ("Normal", "Stage III"),
                                       ("Normal", "Stage IV")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('endometrial cancer')

en_trans_cli_plot(gene = "EIF4E")

## Merge clinical attribute with transcriptomics dataframe
def en_trans_cli_plot(gene):
    en_clinical_and_proteomics = en.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "transcriptomics",
            metadata_cols    = "Proteomics_Tumor_Normal",
            omics_genes      = gene)
    en_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    en_clinical_and_proteomics["Proteomics_Tumor_Normal"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    ax = sns.boxplot(x          = "Proteomics_Tumor_Normal",
                     y          = gene + '_transcriptomics',
                     data       = en_clinical_and_proteomics,
                     showfliers = False)
    sns.stripplot(x        = "Proteomics_Tumor_Normal",
                  y        = gene + '_transcriptomics',
                  data     = en_clinical_and_proteomics,
                  color    = '.3')
    add_stat_annotation(ax,
                        data        = en_clinical_and_proteomics,
                        x           = "Proteomics_Tumor_Normal",
                        y           = gene + '_transcriptomics',
                        boxPairList = [("Tumor", "Adjacent_normal")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('endometrial cancer')

en_trans_cli_plot(gene = "EIF4E")


###################################################################
## Associating Clinical Variables with proteomics in colon cancer##
###################################################################
### load the dataframe for clinical results by calling the en.get_clinical() method
col_clinical_data = col.get_clinical()
print(col_clinical_data.columns)

## Choose Clinical Attribute and Merge Dataframes
col_clinical_attribute = "Stage"

## Merge clinical attribute with proteomics dataframe
## Merge clinical attribute with proteomics dataframe
def col_pro_cli_plot(gene):
    col_clinical_and_proteomics = col.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "proteomics",
            metadata_cols    = "Stage",
            omics_genes      = gene)
    col_clinical_and_proteomics["Stage"] = col_clinical_and_proteomics["Stage"].fillna("Normal")
    col_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    col_clinical_and_proteomics["Stage"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    order      = ["Normal", "Stage I", "Stage II", "Stage III", "Stage IV"]
    ax = sns.boxplot(x          = "Stage",
                     y          = gene + '_proteomics',
                     data       = col_clinical_and_proteomics,
                     showfliers = False,
                     order      = order)
    ax = sns.stripplot(x        = "Stage",
                  y        = gene + '_proteomics',
                  data     = col_clinical_and_proteomics,
                  color    = '.3',
                  order     = order)
    add_stat_annotation(ax,
                        data        = col_clinical_and_proteomics,
                        x           = "Stage",
                        y           = gene + '_proteomics',
                        order       = order,
                        boxPairList = [("Normal", "Stage I"),
                                       ("Normal", "Stage II"),
                                       ("Normal", "Stage III"),
                                       ("Normal", "Stage IV")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('colon cancer')

col_pro_cli_plot(gene = "EIF4E")

###########################################################################
## Associating Clinical Variables with Phosphoproteomics in colon cancer ##
###########################################################################
def col_pho_cliplot(gene):
    col_clinical_and_proteomics = col.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics",
            metadata_cols    = "Stage")
    col_clinical_and_proteomics["Stage"] = col_clinical_and_proteomics["Stage"].fillna("Normal")
    col_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    col_clinical_and_proteomics["Stage"].unique()
    PhosphoSite = list(col_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in PhosphoSite:
        try:
            print(i)
            col_clinical_and_proteomics = col.join_metadata_to_omics(
                metadata_df_name = "clinical",
                omics_df_name    = "phosphoproteomics",
                metadata_cols    = "Stage")
            col_clinical_and_proteomics["Stage"] = col_clinical_and_proteomics["Stage"].fillna("Normal")
            col_clinical_and_proteomics = col_clinical_and_proteomics.dropna(subset = [i])
            plt.figure()
            sns.set(style      ="white",
                    font_scale = 1.5)
            order = ["Normal",
                     "Stage I",
                     "Stage II",
                     "Stage III",
                     "Stage IV"]
            ax = sns.boxplot(x          = "Stage",
                             y          = i,
                             data       = col_clinical_and_proteomics,
                             showfliers = False,
                             order      = order)
            sns.stripplot(x        = "Stage",
                          y        = i,
                          data     = col_clinical_and_proteomics,
                          color    = '.3',
                          order     = order)
            add_stat_annotation(ax,
                                data        = col_clinical_and_proteomics,
                                x           = "Stage",
                                y           = i,
                                order       = order,
                                boxPairList = [("Normal", "Stage I"),
                                               ("Normal", "Stage II"),
                                               ("Normal", "Stage III"),
                                               ("Normal", "Stage IV")],
                                test        = 't-test_ind',
                                textFormat  = 'star',
                                loc         = 'inside',
                                verbose     = 2)
            plt.title('colon cancer')
        except: ValueError
        pass

col_pho_cliplot(gene = "MKNK2_")

#########################################################################
## Associating Clinical Variables with Transcriptomics in colon cancer ##
#########################################################################
def col_tra_cli_plot(gene):
    col_clinical_and_proteomics = col.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "transcriptomics",
            metadata_cols    = "Stage",
            omics_genes      = gene)
    col_clinical_and_proteomics["Stage"] = col_clinical_and_proteomics["Stage"].fillna("Normal")
    col_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    col_clinical_and_proteomics["Stage"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    order = ["Normal",
             "Stage I",
             "Stage II",
             "Stage III",
             "Stage IV"]
    ax = sns.boxplot(x          = "Stage",
                     y          = gene + '_transcriptomics',
                     data       = col_clinical_and_proteomics,
                     showfliers = False,
                     order      = order)
    ax = sns.stripplot(x        = "Stage",
                  y        = gene + '_transcriptomics',
                  data     = col_clinical_and_proteomics,
                  color    = '.3',
                  order     = order)
    add_stat_annotation(ax,
                        data        = col_clinical_and_proteomics,
                        x           = "Stage",
                        y           = gene + '_transcriptomics',
                        order       = order,
                        boxPairList = [("Normal", "Stage I"),
                                       ("Normal", "Stage II"),
                                       ("Normal", "Stage III"),
                                       ("Normal", "Stage IV")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('colon cancer')

col_tra_cli_plot(gene = "MKNK1")

def col_tra_cli_plot(gene):
    col_clinical_and_proteomics = col.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "transcriptomics",
            metadata_cols    = "Stage",
            omics_genes      = gene)
    col_clinical_and_proteomics["Stage"] = col_clinical_and_proteomics["Stage"].fillna("Normal")
    col_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    col_clinical_and_proteomics["Stage"].unique()
    sns.set(style      ="white",
            font_scale = 1.5)
    order = ["Stage I",
             "Stage II",
             "Stage III",
             "Stage IV"]
    ax = sns.boxplot(x          = "Stage",
                     y          = gene + '_transcriptomics',
                     data       = col_clinical_and_proteomics,
                     showfliers = False,
                     order      = order)
    sns.stripplot(x        = "Stage",
                  y        = gene + '_transcriptomics',
                  data     = col_clinical_and_proteomics,
                  color    = '.3',
                  order     = order)
    add_stat_annotation(ax,
                        data        = col_clinical_and_proteomics,
                        x           = "Stage",
                        y           = gene + '_transcriptomics',
                        order       = order,
                        boxPairList = [("Stage I", "Stage II"),
                                       ("Stage I", "Stage III"),
                                       ("Stage I", "Stage IV")],
                        test        = 't-test_ind',
                        textFormat  = 'star',
                        loc         = 'inside',
                        verbose     = 2)
    plt.title('colon cancer')

col_tra_cli_plot(gene = "EIF4E")


#####################################################################
## Associating Clinical Variables with Proteomics in ovarian cancer##
#####################################################################
## load the dataframe for clinical results by calling the en.get_clinical() method
ov_clinical_data = ov.get_clinical()
print(ov_clinical_data.columns)

## Choose Clinical Attribute and Merge Dataframes
ov_clinical_attribute = "Sample_Tumor_Normal"

## Merge clinical attribute with proteomics dataframe
def ov_pro_cli_plot(gene):
    ov_clinical_and_proteomics = ov.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "proteomics",
            metadata_cols    = "Sample_Tumor_Normal",
            omics_genes      = gene)
    ov_clinical_and_proteomics.head()
    cols = []
    count = 1
    for column in ov_clinical_and_proteomics.columns:
        if column == gene + '_proteomics':
           cols.append(gene + '_proteomics'+ str(count))
           count+=1
           continue
        cols.append(column)
    ov_clinical_and_proteomics.columns = cols
## Show possible variations of Histologic_type
    ov_clinical_and_proteomics["Sample_Tumor_Normal"].unique()
    Genes = list(ov_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in Genes:
            print(i)
            ov_clinical_and_proteomics = ov.join_metadata_to_omics(
                        metadata_df_name = "clinical",
                        omics_df_name    = "proteomics",
                        metadata_cols    = "Sample_Tumor_Normal",
                        omics_genes      = gene)
            cols = []
            count = 1
            for column in ov_clinical_and_proteomics.columns:
                if column == gene + '_proteomics':
                   cols.append(gene + '_proteomics'+ str(count))
                   count+=1
                   continue
                cols.append(column)
            ov_clinical_and_proteomics.columns = cols

            ov_clinical_and_proteomics = ov_clinical_and_proteomics.dropna(subset = [i])
            plt.figure()
            sns.set(style      ="white",
                    font_scale = 1.5)
            order      = ["Normal", "Tumor"]
            ax = sns.boxplot(x          = "Sample_Tumor_Normal",
                             y          = i,
                             data       = ov_clinical_and_proteomics,
                             showfliers = False,
                             order      = order)
            ax = sns.stripplot(x        = "Sample_Tumor_Normal",
                          y        = i,
                          data     = ov_clinical_and_proteomics,
                          color    = '.3',
                          order     = order)
            add_stat_annotation(ax,
                                data        = ov_clinical_and_proteomics,
                                x           = "Sample_Tumor_Normal",
                                y           = i,
                                order       = order,
                                boxPairList = [("Normal", "Tumor")],
                                test        = 't-test_ind',
                                textFormat  = 'star',
                                loc         = 'inside',
                                verbose     = 2)
            plt.title('ovarian cancer')


ov_pro_cli_plot(gene = "EIF4E")

############################################################################
## Associating Clinical Variables with Phosphoproteomics in ovarian cancer##
############################################################################
def ov_pho_cli_plot(gene):
    ov_clinical_and_proteomics = ov.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics",
            metadata_cols    = "Sample_Tumor_Normal",
            omics_genes      = gene)
    ov_clinical_and_proteomics.head()
    ov_clinical_and_proteomics = ov_clinical_and_proteomics.loc[:, ~ov_clinical_and_proteomics.columns.duplicated()]
## Show possible variations of Histologic_type
    ov_clinical_and_proteomics["Sample_Tumor_Normal"].unique()
    Genes = list(ov_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in Genes:
            print(i)
            ov_clinical_and_proteomics = ov.join_metadata_to_omics(
                        metadata_df_name = "clinical",
                        omics_df_name    = "phosphoproteomics",
                        metadata_cols    = "Sample_Tumor_Normal",
                        omics_genes      = gene)
            ov_clinical_and_proteomics = ov_clinical_and_proteomics.loc[:, ~ov_clinical_and_proteomics.columns.duplicated()]
            ov_clinical_and_proteomics = ov_clinical_and_proteomics.dropna(subset = [i])
            plt.figure()
            sns.set_style("white")
            order      = ["Normal", "Tumor"]
            ax = sns.boxplot(x          = "Sample_Tumor_Normal",
                             y          = i,
                             data       = ov_clinical_and_proteomics,
                             showfliers = False,
                             order      = order)
            sns.stripplot(x        = "Sample_Tumor_Normal",
                          y        = i,
                          data     = ov_clinical_and_proteomics,
                          color    = '.3',
                          order     = order)
            add_stat_annotation(ax,
                                data        = ov_clinical_and_proteomics,
                                x           = "Sample_Tumor_Normal",
                                y           = i,
                                order       = order,
                                boxPairList = [("Normal", "Tumor")],
                                test        = 't-test_ind',
                                textFormat  = 'star',
                                loc         = 'inside',
                                verbose     = 2)
            plt.title('ovarian cancer')


ov_pho_cli_plot(gene = "EIF4EBP1")



def ovcliplot(gene):
    ov_clinical_and_proteomics = ov.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics",
         #   metadata_cols    = "Tumor_Stage_Ovary_FIGO",
            omics_genes      = gene)
    ov_clinical_and_proteomics["Tumor_Stage_Ovary_FIGO"] = ov_clinical_and_proteomics["Tumor_Stage_Ovary_FIGO"].fillna("Normal")
    ov_clinical_and_proteomics.head()
## Show possible variations of Histologic_type
    ov_clinical_and_proteomics["Tumor_Stage_Ovary_FIGO"].unique()
    PhosphoSite = list(ov_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in PhosphoSite:
        print(i)
        ov_clinical_and_proteomics = ov.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "phosphoproteomics",
         #   metadata_cols    = "Tumor_Stage_Ovary_FIGO",
            omics_genes      = gene)
        ov_clinical_and_proteomics["Tumor_Stage_Ovary_FIGO"] =      ov_clinical_and_proteomics["Tumor_Stage_Ovary_FIGO"].fillna("Normal")
        # ov_clinical_and_proteomics = ov_clinical_and_proteomics.dropna(subset = [i])
        plt.figure()
        sns.set_style("white")
        order = ["Normal",
                 "IIIA",
                 "IIIB",
                 "IIIC",
                 "IV"]
        ax = sns.boxplot(x          = "Tumor_Stage_Ovary_FIGO",
                         y          = i,
                         data       = ov_clinical_and_proteomics,
                         showfliers = False,
                         order      = order)
        sns.stripplot(x        = "Tumor_Stage_Ovary_FIGO",
                      y        = i,
                      data     = ov_clinical_and_proteomics,
                      color    = '.3',
                      order     = order)
        add_stat_annotation(ax,
                            data        = ov_clinical_and_proteomics,
                            x           = "Tumor_Stage_Ovary_FIGO",
                            y           = i,
                            order       = order,
                            boxPairList = [("Normal", "IIIA"),
                                           ("Normal", "IIIB"),
                                           ("Normal", "IIIC"),
                                           ("Normal", "IV")],
                            test        = 't-test_ind',
                            textFormat  = 'star',
                            loc         = 'inside',
                            verbose     = 2)
        plt.title('ovarian cancer')

ovcliplot(gene = "EIF4E")


###########################################################################
## Associating Clinical Variables with Transcriptomics in ovarian cancer ##
###########################################################################
## Merge clinical attribute with proteomics dataframe
def ov_tra_cli_plot(gene):
    ov_clinical_and_proteomics = ov.join_metadata_to_omics(
            metadata_df_name = "clinical",
            omics_df_name    = "transcriptomics",
            metadata_cols    = "Sample_Tumor_Normal",
            omics_genes      = gene)
    ov_clinical_and_proteomics.head()
    cols = []
    count = 1
    for column in ov_clinical_and_proteomics.columns:
        if column == gene + '_transcriptomics':
           cols.append(gene + '_transcriptomics'+ str(count))
           count+=1
           continue
        cols.append(column)
    ov_clinical_and_proteomics.columns = cols
## Show possible variations of Histologic_type
    ov_clinical_and_proteomics["Sample_Tumor_Normal"].unique()
    Genes = list(ov_clinical_and_proteomics.filter(like = gene).columns.values.tolist())
    for i in Genes:
            print(i)
            ov_clinical_and_proteomics = ov.join_metadata_to_omics(
                        metadata_df_name = "clinical",
                        omics_df_name    = "transcriptomics",
                        metadata_cols    = "Sample_Tumor_Normal",
                        omics_genes      = gene)
            cols = []
            count = 1
            for column in ov_clinical_and_proteomics.columns:
                if column == gene + '_transcriptomics':
                   cols.append(gene + '_transcriptomics'+ str(count))
                   count+=1
                   continue
                cols.append(column)
            ov_clinical_and_proteomics.columns = cols

            ov_clinical_and_proteomics = ov_clinical_and_proteomics.dropna(subset = [i])
            plt.figure()
            sns.set(style      ="white",
                    font_scale = 1.5)
            order = ["Normal", "Tumor"]
            ax = sns.boxplot(x          = "Sample_Tumor_Normal",
                             y          = i,
                             data       = ov_clinical_and_proteomics,
                             showfliers = False,
                             order      = order)
            ax = sns.stripplot(x        = "Sample_Tumor_Normal",
                          y        = i,
                          data     = ov_clinical_and_proteomics,
                          color    = '.3',
                          order     = order)
            add_stat_annotation(ax,
                                data        = ov_clinical_and_proteomics,
                                x           = "Sample_Tumor_Normal",
                                y           = i,
                                order       = order,
                                boxPairList = [("Normal", "Tumor")],
                                test        = 't-test_ind',
                                textFormat  = 'star',
                                loc         = 'inside',
                                verbose     = 2)
            plt.title('ovarian cancer')


ov_tra_cli_plot(gene = "EIF4E")












### Extra controls

## correlation between JUN-S243 and SOX2 mRNA in endometrial cancer
def encorrelationplot(gene1, gene2):
    gene_cross_en = en.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_en.head())
    gene_cross_en = gene_cross_en.dropna(subset = ['JUN-S243_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['JUN-S243_phosphoproteomics'],
                    gene_cross_en['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style      ="white",
            font_scale = 1.5)
    plot = sns.regplot(x    = gene_cross_en['JUN-S243_phosphoproteomics'],
                       y    = gene_cross_en['SOX2_transcriptomics'],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-S243 phosphoproteomics',
             ylabel = 'SOX2 transcriptomics',
             title  = 'Proteomics vs Transcriptomics (Endometrial Cancer)')
    plt.show()

encorrelationplot(gene1 = "JUN", gene2 = "SOX2")

## correlation between JUN-T239 and SOX2 mRNA in endometrial cancer
def encorrelationplot(gene1, gene2):
    gene_cross_en = en.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_en.head())
    gene_cross_en = gene_cross_en.dropna(subset = ['JUN-T239_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['JUN-T239_phosphoproteomics'],
                    gene_cross_en['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style      ="white",
            font_scale = 1.5)
    plot = sns.regplot(x    = gene_cross_en['JUN-T239_phosphoproteomics'],
                       y    = gene_cross_en['SOX2_transcriptomics'],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-T239 phosphoproteomics',
             ylabel = 'SOX2 transcriptomics',
             title  = 'Proteomics vs Transcriptomics (Endometrial Cancer)')
    plt.show()

encorrelationplot(gene1 = "JUN", gene2 = "SOX2")

## correlation between JUN-S58 and SOX2 mRNA in endometrial cancer
def encorrelationplot(gene1, gene2):
    gene_cross_en = en.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_en.head())
    gene_cross_en = gene_cross_en.dropna(subset = ['JUN-S58_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['JUN-S58_phosphoproteomics'],
                    gene_cross_en['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style = "white")
    plot = sns.regplot(x    = gene_cross_en['JUN-S58_phosphoproteomics'],
                       y    = gene_cross_en['SOX2_transcriptomics'],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-S58 phosphoproteomics',
             ylabel = gene2 + ' transcriptomics',
             title  = 'Proteomics vs Transcriptomics (Endometrial Cancer)')
    plt.show()

encorrelationplot(gene1 = "JUN", gene2 = "SOX2")

## correlation between JUN-S63 and SOX2 mRNA in endometrial cancer
def encorrelationplot(gene1, gene2):
    gene_cross_en = en.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_en.head())
    gene_cross_en = gene_cross_en.dropna(subset = ['JUN-S63_phosphoproteomics'])
    print(gene_cross_en.head())
    corr = pearsonr(gene_cross_en['JUN-S63_phosphoproteomics'],
                    gene_cross_en['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style = "white")
    plot = sns.regplot(x    = gene_cross_en['JUN-S63_phosphoproteomics'],
                       y    = gene_cross_en['SOX2_transcriptomics'],
                       data = gene_cross_en)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-S63 phosphoproteomics',
             ylabel = gene2 + ' transcriptomics',
             title  = 'Proteomics vs Transcriptomics for ' + gene1 + ' (Endometrial Cancer)')
    plt.show()

encorrelationplot(gene1 = "JUN", gene2 = "SOX2")

## correlation between JUN-S63 and SOX2 mRNA in ovarian cancer
def ovcorrelationplot(gene1, gene2):
    gene_cross_ov = ov.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_ov.head())
    gene_cross_ov = gene_cross_ov.dropna(subset = ['JUN-S63s_phosphoproteomics'])
    print(gene_cross_ov.head())
    corr = pearsonr(gene_cross_ov['JUN-S63s_phosphoproteomics'],
                    gene_cross_ov['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style = "white")
    plot = sns.regplot(x    = gene_cross_ov['JUN-S63s_phosphoproteomics'],
                       y    = gene_cross_ov['SOX2_transcriptomics'],
                       data = gene_cross_ov)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-S63 phosphoproteomics',
             ylabel = gene2 + ' transcriptomics',
             title  = 'Proteomics vs Transcriptomics for ' + gene1 + ' (Ovarian Cancer)')
    plt.show()

ovcorrelationplot(gene1 = "JUN", gene2 = "SOX2")

## correlation between JUN-S73 and SOX2 mRNA in ovarian cancer
def ovcorrelationplot(gene1, gene2):
    gene_cross_ov = ov.join_omics_to_omics(df1_name = "phosphoproteomics",
                                           df2_name = "transcriptomics",
                                           genes1   = gene1,
                                           genes2   = gene2)
    print(gene_cross_ov.head())
    gene_cross_ov = gene_cross_ov.dropna(subset = ['JUN-S73s_phosphoproteomics'])
    #gene_cross_ov = np.log2(gene_cross_ov)
    print(gene_cross_ov.head())
    corr = pearsonr(gene_cross_ov['JUN-S73s_phosphoproteomics'],
                    gene_cross_ov['SOX2_transcriptomics'])
    corr = [np.round(c, 2) for c in corr]
    print(corr)
    sns.set(style = "white")
    plot = sns.regplot(x    = gene_cross_ov['JUN-S73s_phosphoproteomics'],
                       y    = gene_cross_ov['SOX2_transcriptomics'],
                       data = gene_cross_ov)
    text = 'r=%s, p=%s' % (corr[0], corr[1])
    tl   = ((plot.get_xlim()[1] - plot.get_xlim()[0])*0.010 + plot.get_xlim()[0],
            (plot.get_ylim()[1] - plot.get_ylim()[0])*0.95 + plot.get_ylim()[0])
    plot.text(tl[0], tl[1], text, fontsize=12)
    plot.set(xlabel = 'JUN-S73 phosphoproteomics',
             ylabel = gene2 + ' transcriptomics',
             title  = 'Proteomics vs Transcriptomics for ' + gene1 + ' (Ovarian Cancer)')
    plt.show()

ovcorrelationplot(gene1 = "JUN", gene2 = "SOX2")
