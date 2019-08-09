import xenaPython as xena

GENES = ['FOXM1', 'TP53']

def get_codes(host, dataset, fields, data):
    "get codes for enumerations"
    codes = xena.field_codes(host, dataset, fields)
    codes_idx = dict([(x['name'], x['code'].split('\t')) for x in codes if x['code'] is not None])
    for i in range(len(fields)):
        if fields[i] in codes_idx:
            data[i] = [None if v == 'NaN' else codes_idx[fields[i]][int(v)] for v in data[i]]
    return data

def get_fields(host, dataset, samples, fields):
    "get field values"
    data = xena.dataset_fetch(host, dataset, samples, fields)
    return data

def get_fields_and_codes(host, dataset, samples, fields):
    "get fields and resolve codes"
    return get_codes( host, dataset, fields, get_fields( host, dataset, samples, fields))


cohort = 'TCGA TARGET GTEx'
host = xena.PUBLIC_HUBS['toilHub']
# get samples IDs in cohort
samples = xena.cohort_samples(host, cohort, None)
samples[0: 10]
# get expression for GENES
dataset = 'TCGA-GTEx-TARGET-gene-exp-counts.deseq2-normalized.log2'
expression = xena.dataset_gene_probes_values(host, dataset, samples, GENES)
expression = get_fields(host, dataset, samples, GENES) # list of lists.
expression_by_gene = dict([(g['gene'], g['scores'][0]) for g in expression])
[expression_by_gene.keys(), GENES[0], expression_by_gene[GENES[0]][0:10]]


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