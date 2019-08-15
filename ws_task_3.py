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
import re

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
	
	clinvar_variant_dict = {}
	dbsnp_variant_dict = {}
	var_to_remove_list = []
	#put it into a dict

	for dict_1_key, dict_1_value in gene_details.items():
		for dict_2_key, dict_2_value in dict_1_value.items():
			# dict 2 value is a list of dicts containing individual entries
			if type(dict_2_value) == list:
				# capture different variables from clinvar or dbsnp
				for dict_3 in dict_2_value:
					# dict3 contains mutiple dicts
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

	#create list of variants to remove from the dictionary
	for variant_key,info_value in clinvar_variant_dict.items():
		if info_value['rs_number'] == 'rsNA' or 'NM_' not in variant_key:
			var_to_remove_list.append(variant_key)
	# remove variants from dictionary
	for variant in var_to_remove_list:
		clinvar_variant_dict.pop(variant)
	return clinvar_variant_dict


gene_dict = {}


for gene in string_list_of_genes:
	term = gene + '[gene]'
	clinvar_dict = get_variant_info('clinvar', term, '100000')
	gene_dict.update({gene:clinvar_dict} )	

# writing all to .txt file
for gene, variant_value in gene_dict.items():
	with open('{}_snp_list.txt'.format(gene), 'w') as file:
		file.write("{}\t{}\t{}\t{}\n".format(
			'Variant',
			'Variant type',
			'rs_number',
			'Pathogenicity'
			)
		)
		for variant_key,variant_details in variant_value.items():
				file.write("{}\t{}\t{}\t{}\n".format(
				variant_key,
				variant_details['variant_type'],
				variant_details['rs_number'],
				variant_details['pathogenicity']
				)	
			)