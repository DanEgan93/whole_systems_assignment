'''
	Whole Systems Assignment

	Task 2- Write a script to obtain a list of genes that interact with the
	SUFU gene

'''

from Bio import Entrez
import csv

Entrez.email = "dan_egan@live.co.uk"

handle = Entrez.esearch(
	db = 'gene', 
	term = 'SUFU[ALL Fields] AND (Homo Sapiens [Organism] OR homo\
	 sapiens[ALL Fields]) AND alive[prop]',
	retmax = '100000'
	)

data = Entrez.read(handle)
handle.close()

sufu_list = data['IdList']
sufu_list_query = ",".join(map(str,sufu_list))

handle = Entrez.esummary(db='gene', id=sufu_list_query, retmode='xml')
sufu_details = Entrez.read(handle)

entrez_list_of_genes = []

# Returns in Dict(Dict(list_of_dict(dict)))

for k,v in sufu_details.items():
	for e,i in v.items():
		if type(i) == list:
			for g in i:
				entrez_list_of_genes.append(g['Name'].encode('UTF-8'))

#comparing list of genes from entrez and string
string_list_of_genes = []

with open('string_interactions.tsv', 'r') as file:
	tsvin = csv.reader(file, delimiter='\t')
	for row in tsvin:
		if row[0] == '#node1':
			continue
		string_list_of_genes.append(row[0])
		string_list_of_genes.append(row[1])

string_list_of_genes = set(string_list_of_genes)
entrez_list_of_genes = set(entrez_list_of_genes)

# Compating lists of genes
gene_list_overlap = []
gene_list_non_overlap = []

print(len(entrez_list_of_genes))
print(len(string_list_of_genes))

#Adding genes in list to text file
with open('entrez_gene_list.txt', 'w') as file:
	for gene in entrez_list_of_genes:
		file.write(gene+'\n')
