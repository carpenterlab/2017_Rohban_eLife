# Use this to create yml file

import csv
import sys

#fname = '../data/TA-POP-B1/mean-all_controls-pop_genes.csv'
fname = '../data/TA-POP-B1/mean-all_controls-all_genes.csv'
CTRL = "CTRL"

gene_shrna = dict()
shrna_gene = dict()
gene_shrnaid = dict()
shrna_id = dict()
shrna_role = dict()
id_shrna = dict()

f = open(fname, 'rU')
try:
    reader =  csv.DictReader(f)
    for row in reader:
        gene_shrna.setdefault(row['Gene'], set()).add(row['RNAi'])
        shrna_gene[row['RNAi']] = row['Gene']
        shrna_role[row['RNAi']] = row['Role']
finally:
    f.close()

for shrna, gene in shrna_gene.items():
    id = gene_shrnaid.setdefault(gene, 0)
    id = id + 1
    id_ = '%s.%d' % (gene, id)
    shrna_id[shrna] = id_
    id_shrna[id_] = shrna
    gene_shrnaid[gene] = id

idl = shrna_id.values()
idl.sort()

treatment_names = []
treatment_names_abbrev = []
treatment_names_grouped = []
control_names = []
control_names_abbrev = []
control_names_grouped = []

for id in idl:
    shrna = id_shrna[id]
    if (shrna_role[shrna] != CTRL):
        treatment_names.append(shrna)
        treatment_names_abbrev.append(shrna_id[shrna])
        treatment_names_grouped.append(shrna_gene[shrna])
    else:
        control_names.append(shrna)
        control_names_abbrev.append(shrna_id[shrna])
        control_names_grouped.append(shrna_gene[shrna])

print "treatment_names : [" + ",".join(treatment_names) + "]"
print "treatment_names_abbrev : [" + ",".join(treatment_names_abbrev) + "]"
print "treatment_names_grouped_abbrev : [" + ",".join(treatment_names_grouped) + "]"
print "control_name : [" + ",".join(control_names) + "]"
# Hack ! control_names_abbrev is technically the right thing to place 
# here but the R code requires it to be grouped. 
print "control_name_abbrev : [" + ",".join(control_names_grouped) + "]" 


    
    


