import xenaPython as xena
import pandas as pd

GENES = ['FOXM1', 'TP53']

def get_codes(host, dataset, fields, data):
    "get codes for enumerations"
    codes = xena.field_codes(host, dataset, fields)
    codes_idx = dict([(x['name'], 
                       x['code'].split('\t')) for x in codes if x['code'] is not None])
    for i in range(len(fields)):
        if fields[i] in codes_idx:
            data[i] = [None if v == 'NaN' else codes_idx[fields[i]][int(v)] for v in data[i]]
    return data


def get_fields(host, dataset, samples, fields):
    "get field values, column names in the spreadsheet"
    data = xena.dataset_fetch(host, dataset, samples, fields)
    return data


def get_fields_and_codes(host, dataset, samples, fields):
    "get fields and resolve NA in the value"
    return get_codes(host, dataset, fields, get_fields( host, dataset, samples, fields))

# dictionary with all hub links
xena.PUBLIC_HUBS  
# pancanAtlas cohort
cohort = 'TCGA PanCanAtlas'
host = xena.PUBLIC_HUBS['pancanAtlasHub']

    
# get expression for GENES
expression_dataset = 'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena'
samples = xena.dataset_samples(host, expression_dataset, None)
samples[0: 10]
expression = get_fields_and_codes(host, 
                                  expression_dataset, 
                                  samples, 
                                  GENES) # list of lists.
expression_by_gene = dict(zip(GENES, expression))      # index by gene.
[expression_by_gene.keys(), GENES[0], expression_by_gene[GENES[0]][0:10]]
# note that missing data is returned as 'NaN'. One might want to remap this to None or NaN, depending on the later analysis tools.


# get disease type and survival columns
survival_dataset = 'Survival_SupplementalTable_S1_20171025_xena_sp'
fields = ['cancer type abbreviation', 'OS', 'OS.time']
values = get_fields_and_codes(host, 
                              survival_dataset, 
                              samples, 
                              fields) # list of lists
phenotypes = dict(zip(fields, values)) # index by phenotype
# show all unique variable in the list phenotypes['cancer type abbreviation']
phenotype_index = set(phenotypes['cancer type abbreviation'])
print(phenotype_index)


# get sample type. TCGA includes a few "normal" tissue samples. These normals are of
# limited value because there are few of them, and they are not entirely normal, being
# taken from disease tissue, outside of the visible tumor. It's often best to omit them.
sampletype_dataset = 'TCGA_phenotype_denseDataOnlyDownload.tsv'
fields = ['sample_type']
values = get_fields_and_codes(host, 
                              sampletype_dataset, 
                              samples, 
                              fields)
set(values[0])



## has to use ['None'] to list all cohorts within the hub
## all_cohorts(host, exclude)
xena.all_cohorts(hub, ['None']) 
cohort = 'https://pancanatlas.xenahubs.net'

xena.cohort_summary(hub, ['None']) # count datasets per cohort
## Dataset metadata for datasets in the given cohorts
xena.dataset_list(hub, ['TCGA Pan-Cancer (PANCAN)'])
## use dataset ID in Xena website
dataset = 'EB++AdjustPANCAN_IlluminaHiSeq_RNASeqV2.geneExp.xena'
gene = ["VEGFA"]
## samples = xena.dataset_samples(hub, dataset, None)
samples = xena.dataset_samples(hub, dataset, 10)
## get expression for gene
expression = xena.dataset_gene_probes_values(hub, dataset, samples, gene)
expression_by_gene = dict([(g['gene'], g['scores'][0]) for g in expression])
[expression_by_gene.keys(), gene[0], expression_by_gene[gene[0]][0:10]]

## retrieve sample type value
## NOTE: different hub
hub = 'https://pancanatlas.xenahubs.net'
dataset = 'TCGA_phenotype_denseDataOnlyDownload.tsv'
fields = ['_primary_disease', 'sample_type']
# _sample_type will identify normals. 
# _primary_disease will identify cancer study group
values = get_fields_and_codes(hub, dataset, samples, fields) # list of lists
phenotypes = dict(zip(fields, values)) # index by phenotype
phenotypes['_primary_disease'][0:10]



toil_summary = {
    'samples': samples,
    'expression': expression_by_gene,
    'phenotypes': phenotypes
}




# list of field (probes) and probeMap, 
# a probeMap is used to map a gene location to a set
# of probes (all RNA-seq counts at a given genomic location)
xena.dataset_field(hub, dataset)

# Exon counts in gene, and probe genomic positions, for given samples
lst = xena.dataset_gene_probes_values(hub, dataset, samples, gene)
# create a Pandas dataframe from probes_value list
data1 =  pd.DataFrame(lst[0])
data2 =  pd.DataFrame(lst[1])
data3 = pd.concat([data1, data2], axis=1)
# iris3 = iris3.T


[position, [VEGFA]] = xena.dataset_probe_values(hub, dataset, samples, gene)

# use metadata function to find annotation file for genes covering all field (probes)
xena.dataset_metadata(hub, dataset)
# the probeMap for exon data is 'unc_v2_exon_hg19_probe_TCGA'
probemap = 'unc_v2_exon_hg19_probe_TCGA'







dataset = 'TcgaTargetGTEX_phenotype.txt'
fields = ['_study', '_sample_type']
# As in pancan, there are normal samples in tcga which should probably be removed. _sample_type will
# identify normals. _study will identify tcga vs. gtex vs. target.
values = get_fields_and_codes(host, dataset, samples, fields) # list of lists
phenotypes = dict(zip(fields, values)) # index by phenotype
phenotypes['_study'][0:10]


toil_summary = {
    'samples': samples,
    'expression': expression_by_gene,
    'phenotypes': phenotypes
}