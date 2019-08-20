'''
	Whole Systems Assignment

	Task 4- Write a script to:
	1) Obtain OMIM code for Meier-Gorlin Syndrome 1 and write this to a text file. 
	2) Create a list of papers linked to the OMIM code of the diesease/condition
	of interest and write this to a text file
	3) Create a list of all papers that mention Meier-Gorlin Syndrome and write this 
	to a text file.
	4) Crate a list of all SNPs in ClinVar related Meir-Gorlin Syndrome and write this
	to a text file.
	
	Note- For this assignment I have chosen the MEIER-GORLIN SYNDROME 1.
	The reason for this is that there are many types of the sydrome with
	each displaying slightly different phenotypes and type 1 is the most
	common version of the syndrome.

'''

from Bio import Entrez

Entrez.email = 'dan_egan@live.co.uk'



def get_omim_code(database_name, term, retmax, disease_search_term):
	'''
	The get_omim_code functions uses Entrez's esearch method to returns the
	number id of omim entries linked to the term 'Meier-Gorlin Syndrome' and 
	collects the omim code for the MEIER-GORLIN SYNDROME 1; MGORS1 term only.
	See above for description. 

	'''

	handle = Entrez.esearch(db =database_name, term=term, retmax=retmax)
	data = Entrez.read(handle)

	omim_gene_list = data['IdList']
	omim_gene_list = ",".join(map(str,omim_gene_list))

	handle = Entrez.esummary(db='OMIM', id=omim_gene_list, retmax='100000')
	omim_details = Entrez.read(handle)

	for omim_entry in omim_details:
		if omim_entry['Title'] == disease_search_term:
			omim_code = omim_entry['Oid'].strip('#')
	return omim_code

def get_omim_linked_papers(db, dbfrom, omim_id):
	'''
	This function uses Entrez's elink method to return the pubmed IDs of
	papers associated with the omim code passed to this function. The esummary
	method is then used to retrieve paper details from the IDs provided.
	'''
	handle = Entrez.elink(db=db,dbfrom=dbfrom, id=omim_id)
	record = Entrez.read(handle)

	papers_linked_to_omim_id = []

	for link in record[0]['LinkSetDb'][-1]['Link']:
		papers_linked_to_omim_id.append(link['Id'])

	papers_linked_to_omim_id = ",".join(map(str,papers_linked_to_omim_id))

	handle = Entrez.esearch(db='pubmed', term=papers_linked_to_omim_id, retmax='100000')
	data = Entrez.read(handle)
	pubmed_search_id_list = data['IdList']

	pubmed_search_id_list = ",".join(map(str,pubmed_search_id_list))
	handle = Entrez.esummary(db='pubmed', id=pubmed_search_id_list, retmax='100000')
	test_deets = Entrez.read(handle)

	gorlin_paper_details = {}

	for record in test_deets:
		# create comma separated list of authors
		author_list = [] 
		for name in record['AuthorList']:
			author_list.append(name.encode('UTF-8'))
		author_string = ', '.join(map(str,author_list)).rstrip(",")
		title = str(record['Title'].encode('UTF-8'))
		year = str(record['PubDate'][0:4].encode('UTF-8'))
		journal = str(record['FullJournalName'].encode('UTF-8'))
		pages = str(record['Pages'].encode('UTF-8'))

		gorlin_paper_details.update({record['Id']: {
			'Title': title, 
			'Pulication_date': year,
			'Author': author_string,
			'Journal': journal,
			'Pages': pages
			}}
		)

	return gorlin_paper_details

def get_generic_paper_list(db, term, retmax):
	''' The esearch and esummary methods are used to return a list of papers
	and paper detailed that mention the condition Meir-Gorlin syndrome 
	(non-syndrome number specific) '''

	handle = Entrez.esearch(db =db, term=term, retmax=retmax)
	data = Entrez.read(handle)

	gorlin_gen_papers = data['IdList']
	gorlin_gen_papers = ",".join(map(str,gorlin_gen_papers))

	handle = Entrez.esummary(db='pubmed', id=gorlin_gen_papers, retmax='100000')
	gorlin_gen_details = Entrez.read(handle)

	gorlin_gen_dict = {}

	for record in gorlin_gen_details:
		# create comma separated list of authors
		author_list = [] 
		for name in record['AuthorList']:
			author_list.append(name.encode('UTF-8'))
		author_string = ', '.join(map(str,author_list)).rstrip(",")
		title = str(record['Title'].encode('UTF-8'))
		year = str(record['PubDate'][0:4].encode('UTF-8'))
		journal = str(record['FullJournalName'].encode('UTF-8'))
		pages = str(record['Pages'].encode('UTF-8'))


		gorlin_gen_dict.update({record['Id']: {
			'Title': title, 
			'Pulication_date': year,
			'Author': author_string,
			'Journal': journal,
			'Pages': pages
			}}
		)
	return gorlin_gen_dict

def get_variant_info(database_name, search_term, retmax):
	'''
	This function does the following:
		1) Get ids from the ClinVar database using the esearch functions related
		to search term 'Meier-Gorlin Syndrome'
		2) Uses esummary to extract variant info 
		3) Add info to a dictionary-  clinvar_variant_dict
		4) Returns clinvar_variant_dict
	'''
	
	handle = Entrez.esearch(
		db = database_name, 
		term = search_term,
		retmax = retmax
		)

	data = Entrez.read(handle)
	handle.close()

	clinvar_entries_id_list = data['IdList']
	clinvar_entries_query = ",".join(map(str,clinvar_entries_id_list))

	handle = Entrez.esummary(db=database_name, id=clinvar_entries_query, retmode='xml')
	clinvar_entry_details = Entrez.read(handle)
	
	clinvar_variant_dict = {}
	dbsnp_variant_dict = {}
	var_to_remove_list = []

	'''
	Esummary function returns gene associated data in the 
	dict_1(dict_2(list_of_dicts(dict_3(dict_4)))) format. The for loop below itterates
	through this structure to extract SNP information for each gene returned from
	the search term passed to the get_variant_info function. 

	'''
	for dict_1_key, dict_1_value in clinvar_entry_details.items():
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


term = 'Meier-Gorlin Syndrome'
disease_search_term = 'MEIER-GORLIN SYNDROME 1; MGORS1'

omim_code = get_omim_code(database_name='omim',term=term, retmax='100000', disease_search_term=disease_search_term)
omim_linked_paper_dict = get_omim_linked_papers(db='pubmed',dbfrom='omim',omim_id=omim_code)
generic_paper_list = get_generic_paper_list(db='pubmed', term=term, retmax='100000')
var_info = get_variant_info(database_name='clinvar', search_term=term, retmax='100000')


'''
The following section writes all information collected in the function above to
independent text files. The names and content of theses files are as follows:
	1) 'omim_code.txt'- OMIM code associated with Meier-Gorlin syndrome 1.
	2) 'omim_paper_list.txt'- List of papers in pubmed linked to above OMIM code.
	3) 'gorlin_generic_search_paper_list.txt'- List of papers in pubmed associated with
	the generic term 'Meier-Gorlin Syndrome'.
	4) 'meier-Gorlin_SNP_list.txt'- List of SNPS associated with Meier-Gorlin Syndrome.
'''

with open('omim_code.txt', 'w') as file:
	file.write('Disease/condition:{}\n'.format(disease_search_term))
	file.write('OMIM code:{}\n'.format(omim_code))

with open("omim_paper_list.txt","w") as file:
	# Creating header
	file.write("{}\t{}\t{}\t{}\t{}\n".format(
				'Title',
				'Pulication Year',
				'Author(s)',
				'Journal',
				'Pages')
			)
	# Adding paper information
	for k, v in omim_linked_paper_dict.items():
		file.write("{}\t{}\t{}\t{}\t{}\n".format(
			v['Title'],
			v['Pulication_date'],
			v['Author'],
			v['Journal'],
			v['Pages'])
		)

with open("gorlin_generic_search_paper_list.txt","w") as file:
	# Creating header
	file.write("{}\t{}\t{}\t{}\t{}\n".format(
				'Title',
				'Pulication Year',
				'Author(s)',
				'Journal',
				'Pages')
			)
	# Adding paper information
	for k, v in generic_paper_list.items():
		file.write("{}\t{}\t{}\t{}\t{}\n".format(
			v['Title'],
			v['Pulication_date'],
			v['Author'],
			v['Journal'],
			v['Pages'])
		)

with open('meier-Gorlin_SNP_list.txt', 'w') as file:
	file.write("{}\t{}\t{}\t{}\n".format(
		'Variant',
		'Clinical significance',
		'Phenotype',
		'rsID',
		'Variant type'
		)
	)
	for variant_key,variant_details in var_info.items():
			file.write("{}\t{}\t{}\t{}\n".format(
			variant_key,
			variant_details['pathogenicity'],
			variant_details['phenotype'],
			variant_details['rs_number'],
			variant_details['variant_type']
			)	
		)