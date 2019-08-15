'''
	Whole Systems Assignment

	Task 2- Using the list of genes from task 2 extract information from ClinVar 
	and dbSNP. For each gene in your list, create a text file containing:
		1) All the SNPs within those genes that are in CLinVar AND dbSNP
		2) The clinical significance of each of the SNPs 
		3) The disease or phenotype with which each of the SNPs is associated

	Note- The list of genes used for this excercise will be the genes extracted
	from the STRING database. This decision was made as this will limit the number
	or requests to the Entrez API.
'''

from Bio import Entrez
import csv

Entrez.email = 'dan_egan@live.co.uk'


''' Creates a list of genes from the string_interactions.tsv file downloaded
String databse.
'''
string_list_of_genes = []

with open('string_interactions.tsv', 'r') as file:
	tsvin = csv.reader(file, delimiter='\t')
	for row in tsvin:
		if row[0] == '#node1':
			continue
		string_list_of_genes.append(row[0])
		string_list_of_genes.append(row[1])
string_list_of_genes = set(string_list_of_genes)


'''
iterate through list of genes
	1) get ids related to search term (clinvar)
	2) uses esummary to exract variant info
	3) repeat for (dbsnp)
	4) combine the lists
'''
	
def get_variant_info(database_name, search_term, retmax):
	
	handle = Entrez.esearch(
		db = database_name, 
		term = search_term,
		retmax = retmax
		)

	data = Entrez.read(handle)
	handle.close()

	gene_id_list = data['IdList']
	gene_id_query = ",".join(map(str,gene_id_list))

	handle = Entrez.esummary(db=database_name, id=gene_id_query, retmode='xml')
	gene_details = Entrez.read(handle)
	
	hgvs_list = []
	pathogenicity_list = []
	pheno_list = []
	variant_type_list = []

	clinvar_variant_dict = {}
	#put it into a dict

	for dict_1_key, dict_1_value in gene_details.items():
		for dict_2_key, dict_2_value in dict_1_value.items():
			if type(dict_2_value) == list:
				for list_of_dict_3 in dict_2_value:
					for dict_4 in list_of_dict_3['variation_set']:
						hgvs = dict_4['cdna_change']
						#need to remove non NM's 
						variant_type = dict_4['variant_type']
					pathogenicity = list_of_dict_3['clinical_significance']['description']
					phenotype = list_of_dict_3['trait_set'][0]['trait_name']

					clinvar_variant_dict.update({hgvs: {
						'variant_type': variant_type, 
						'pathogenicity': pathogenicity,
						'phenotype': phenotype
						}}
					)
					

					# for k,v in list_of_dict_3.items():
					# 	print(k)
						# Variant (HGVS)
						# Clinical_sig
						# Phenotype
						# Variant_type

	return clinvar_variant_dict


clinvar_dict = get_variant_info('clinvar', 'SUFU[gene]', '100000')
print(clinvar_dict)
# get_variant_info('snp', 'SUFU[gene]','100000')