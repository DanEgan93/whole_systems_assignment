'''
	Whole Systems Assignment

	Task 1- Write a script to obtain a list of papers that mention the SUFU gene
	        (including paper details e.g. title, year of publication, author(s),
	        journal of publication and pages). The output from the script must be
	        a .txt file containing the list of papers and paper details.

'''

from Bio import Entrez

Entrez.email = "dan_egan@live.co.uk"

# esearch- returns a list of ids for a pubmed search.
handle = Entrez.esearch(
	db='pubmed',
	retmax= '100000',
	retmode='xml',
	term='SUFU AND gene'
	)
data = Entrez.read(handle)
sufu_list = data['IdList']
handle.close()

# esummary-returns information (inc title) regarding the PubMed IDs provided.
sufu_list_query = ",".join(map(str,sufu_list))

paper_list = []
handle = Entrez.esummary(db='pubmed', id=sufu_list_query, retmode='xml')
sufu_details = Entrez.read(handle)

sufu_paper_details = {}


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

# Write list of publications to gene_list.txt file
with open("paper_list.txt","w") as file:
	# Creating header
	file.write("{}\t{}\t{}\t{}\t{}\n".format(
				'Title',
				'Pulication Year',
				'Author(s)',
				'Journal',
				'Pages')
			)
	# Adding 
	for k, v in sufu_paper_details.items():
		file.write("{}\t{}\t{}\t{}\t{}\n".format(
			v['Title'],
			v['Pulication_date'],
			v['Author'],
			v['Journal'],
			v['Pages'])
		)
file.close()
