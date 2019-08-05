'''
	Whole Systems Assignment

	Task 1- Write a script to obtain a list of papers that mention the SUFU gene
	        (including paper details e.g. title, year of publication, author(s),
	        journal of publication and pages). The output from the script must be
	        a .txt file containing the list of papers and paper details.

'''

from Bio import Entrez

Entrez.email = "dan_egan@live.co.uk"

'''
The esearch function from the Enzrez API has been used used to return the
pubmed IDs (max 1000) of all entries in the database which include the
term 'SUFU AND gene'.

'''
handle = Entrez.esearch(
	db='pubmed',
	retmax= '100000',
	retmode='xml',
	term='SUFU AND gene'
	)
data = Entrez.read(handle)
# Extacting the PubMed ids only
sufu_list = data['IdList']
handle.close()

''' 
Joining the list of PubMed ids into a string (required for the Entrez
esummary function)
'''
sufu_list_query = ",".join(map(str,sufu_list))

paper_list = []

'''
The esummary function will return summary information for all PubMed IDs passed
in the sufu_list_query. This includes paper details e.g. title, 
year of publication, authors... 
'''

handle = Entrez.esummary(db='pubmed', id=sufu_list_query, retmode='xml')
sufu_details = Entrez.read(handle)

sufu_paper_details = {}

'''
Iterating through sufu_details, a list of dictionaries. Each dictionary or 'record' has a
series of information about a single paper. These are accesed e.g. record['title']
and are added to a new dictionary, sufu_paper_details.
'''
for record in sufu_details:
	# create comma separated list of authors
	author_list = [] 
	for name in record['AuthorList']:
		author_list.append(name.encode('UTF-8'))
	author_string = ', '.join(map(str,author_list)).rstrip(",")
	title = str(record['Title'].encode('UTF-8'))
	year = str(record['PubDate'][0:4].encode('UTF-8'))
	journal = str(record['FullJournalName'].encode('UTF-8'))
	pages = str(record['Pages'].encode('UTF-8'))

	sufu_paper_details.update({record['Id']: {
		'Title': title, 
		'Pulication_date': year,
		'Author': author_string,
		'Journal': journal,
		'Pages': pages
		}}
	)

'''
A text file, paper_list.txt, is created and a header is added.
Each item from the sufu_paper_details is then added to the text file.
'''
with open("paper_list.txt","w") as file:
	# Creating header
	file.write("{}\t{}\t{}\t{}\t{}\n".format(
				'Title',
				'Pulication Year',
				'Author(s)',
				'Journal',
				'Pages')
			)
	# Adding paper information
	for k, v in sufu_paper_details.items():
		file.write("{}\t{}\t{}\t{}\t{}\n".format(
			v['Title'],
			v['Pulication_date'],
			v['Author'],
			v['Journal'],
			v['Pages'])
		)
file.close()


