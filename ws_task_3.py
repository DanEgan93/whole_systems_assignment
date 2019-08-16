'''
	Whole Systems Assignment

	Task 3- Using the list of genes from task 2 extract information from ClinVar 
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
STRING database.
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



def get_variant_info(database_name, search_term, retmax):
	'''
	Iterate through list of STRING genes
		1) Get ids from the ClinVar database using the esearch functions related
		to search term '*gene_name*[gene]'
		2) Uses esummary to exract variant info (see desciption above)
		3) Add info to a dictionary- clinvar_variant_dict
		4) Return clinvar_variant_dict
	'''
	
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
	
	clinvar_variant_dict = {}
	dbsnp_variant_dict = {}
	var_to_remove_list = []


	'''
	Esummary function returns gene associated data in the 
	dict_1(dict_2(list_of_dicts(dict_3(dict_4)))) format. The for loop below itterates
	through this structure to extract the gene name for each gene returned from
	the search term described above. 

	'''
	for dict_1_key, dict_1_value in gene_details.items():
		for dict_2_key, list_of_dicts in dict_1_value.items():
			if type(list_of_dicts) == list:
				for dict_3 in list_of_dicts:
					for dict_4 in dict_3['variation_set']:
						if len(dict_4['variation_xrefs']) != 0:
							for source in dict_4['variation_xrefs']:
								if source['db_source'] == 'dbSNP':
									rs_id = source['db_id']
								else:
									rs_id = 'NA'
						else:
							rs_id = 'NA'
						hgvs = dict_4['cdna_change']
						variant_type = dict_4['variant_type']
					pathogenicity = dict_3['clinical_significance']['description']
					phenotype = dict_3['trait_set'][0]['trait_name']
					clinvar_variant_dict.update({hgvs: {
					'variant_type': variant_type, 
					'pathogenicity': pathogenicity,
					'phenotype': phenotype,
					'rs_number': 'rs'+str(rs_id)
					}}
				)

	'''Creates a list of variants to remove from the dictionary including non-
	SNPs and entries that do not have an rs number associated with them i.e.
	the are not present in the dbSNP database.
	'''
	for variant_key,info_value in clinvar_variant_dict.items():
		if info_value['rs_number'] == 'rsNA' or 'NM_' not in variant_key:
			var_to_remove_list.append(variant_key)
	for variant in var_to_remove_list:
		clinvar_variant_dict.pop(variant)
	return clinvar_variant_dict

# Creates gene_dict to store all SNPs from each gene in gene list
gene_dict = {}

# Calls the get_variant_info for each gene in the STRING gene list
for gene in string_list_of_genes:
	term = gene + '[gene]'
	clinvar_dict = get_variant_info('clinvar', term, '100000')
	gene_dict.update({gene:clinvar_dict} )	

# Writing SNP details for each gene to an independent .txt files
for gene, variant_value in gene_dict.items():
	with open('{}_snp_list.txt'.format(gene), 'w') as file:
		file.write("{}\t{}\t{}\t{}\n".format(
			'Variant',
			'Clinical significance',
			'Phenotype',
			'rsID',
			'Variant type'
			)
		)
		for variant_key,variant_details in variant_value.items():
				file.write("{}\t{}\t{}\t{}\n".format(
				variant_key,
				variant_details['pathogenicity'],
				variant_details['phenotype'],
				variant_details['rs_number'],
				variant_details['variant_type']
				)	
			)