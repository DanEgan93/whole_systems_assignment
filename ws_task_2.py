'''
	Whole Systems Assignment

	Task 2- Write a script to obtain a list of genes that interact with the
	SUFU gene

'''

from Bio import Entrez

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

list_of_genes = []

# Returns in Dict(Dict(list_of_dict(dict)))

for k,v in sufu_details.items():
	for e,i in v.items():
		if type(i) == list:
			for g in i:
				list_of_genes.append(g['Name'].encode('UTF-8'))

#comparing list of genes from entrez and string


#Adding genes in list to text file
with open('entrez_gene_list.txt', 'w') as file:
	for gene in list_of_genes:
		file.write(gene+'\n')
