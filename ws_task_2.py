'''
	Whole Systems Assignment

	Task 2- Write a script to obtain a list of genes from the Entrez Gene 
	database which interact with the SUFU gene.

'''

from Bio import Entrez
import csv

Entrez.email = "dan_egan@live.co.uk"

'''
The esearch function from the Enzrez API has been used used to return the
pubmed IDs (max 1000) of all entries in the gene database which include the
term 'SUFU[ALL Fields] AND (Homo Sapiens [Organism] OR 
homo sapiens[ALL Fields]) AND alive[prop]'.

'''

handle = Entrez.esearch(
	db = 'gene', 
	term = 'SUFU[ALL Fields] AND (Homo Sapiens [Organism] OR homo\
	 sapiens[ALL Fields]) AND alive[prop]',
	retmax = '100000'
	)

data = Entrez.read(handle)
handle.close()


''' 
1) Extracting the PubMed ids only
2) Joining the list of PubMed ids into a string (required for the Entrez
esummary function)
'''

sufu_list = data['IdList']
sufu_list_query = ",".join(map(str,sufu_list))

'''
The esummary function will return all genes and associated gene information
using the search term above.
'''
handle = Entrez.esummary(db='gene', id=sufu_list_query, retmode='xml')
sufu_details = Entrez.read(handle)

entrez_list_of_genes = []
string_list_of_genes = []
gene_list_overlap = []
unique_to_string = []
unique_to_entrez = []
genes_to_remove = []

''' Esummary function returns gene associated data in the 
dict_1(dict_2(list_of_dict(dict_3))) format. The for loop below itterates
through this structure to extract the gene name for each gene returned from
the search term described above. 
'''
for dict_1_key, dict_1_value in sufu_details.items():
	for dict_2_key,dict_2_value in dict_1_value.items():
		if type(dict_2_value) == list:
			for dict_3 in dict_2_value:
				entrez_list_of_genes.append(dict_3['Name'].encode('UTF-8'))

''' Creates a list of genes from the string_interactions.tsv file downloaded
String databse.
'''
with open('string_interactions.tsv', 'r') as file:
	tsvin = csv.reader(file, delimiter='\t')
	for row in tsvin:
		if row[0] == '#node1':
			continue
		string_list_of_genes.append(row[0])
		string_list_of_genes.append(row[1])
string_list_of_genes = set(string_list_of_genes)
entrez_list_of_genes = set(entrez_list_of_genes)

'''
Compares the list of genes extracted from the Entrez database and the
list of genes extracted from the String .tsv file provided. A list of genes 
that are common to both lists and two lists containing genes which are unique 
to each data source are then defined. 
'''

for gene in string_list_of_genes:
	if gene in entrez_list_of_genes:
		gene_list_overlap.append(gene)
	else:
		unique_to_string.append(gene)

# Identification and removal of hypothetical genes and miRNA genes.
for gene in entrez_list_of_genes:
	if '_' in gene or 'MIR' in gene:
		genes_to_remove.append(gene)
	gene = gene.strip(' ')
	if gene not in string_list_of_genes:
		unique_to_entrez.append(gene)

for gene in genes_to_remove:
	entrez_list_of_genes.remove(gene)

# The list of genes identified from the Entrez database is written to a file
with open('entrez_gene_list.txt', 'w') as file:
	for gene in entrez_list_of_genes:
		file.write(gene+'\n')

# A text file summarsing the lists of genes from both databases is created.
with open('gene_list_comparison.txt', 'w') as file:
	file.write('*** Total number of genes comparison *** \n')
	file.write('Entrez = {} \n'.format(len(entrez_list_of_genes)))
	file.write('String = {} \n'.format(len(string_list_of_genes)))
	file.write('Number of overlapping genes = {} \n'.format(len(gene_list_overlap)))
	file.write('\n\n')
	file.write('Entrez gene list: \n')
	for gene in entrez_list_of_genes:
		file.write(gene+'\t')
	file.write('\n\n')
	file.write('String gene list: \n')
	for gene in string_list_of_genes:
		file.write(gene+'\t')
	file.write('\n\n')
	file.write('Overlapping genes from both lists:\n')
	for gene in gene_list_overlap:
		file.write(gene+'\t')