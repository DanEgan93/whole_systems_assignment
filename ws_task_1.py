'''
	Whole Systems Assignment

	Task 1- Write a script to obtain a list of papers including the SUFU gene
    		from the PubMed.The output from the script must be a .txt file 
    		containing the list of papers.

'''

from Bio import Entrez

Entrez.email = "dan_egan@live.co.uk"

# esearch- returns a list of ids for a pubmed search.
handle = Entrez.esearch(db='pubmed', retmax= '100000', retmode='xml', term='SUFU')
data = Entrez.read(handle)
sufu_list = data['IdList']
handle.close()

# esummary-returns information (inc title) regarding the PubMed IDs provided.
sufu_list = ",".join(map(str,sufu_list))
paper_list = []
handle = Entrez.esummary(db='pubmed', id=sufu_list, retmode='xml')
sufu_details = Entrez.read(handle)
for sub_record in sufu_details:
	paper_list.append(sub_record['Title'].__repr__())
handle.close()

# Write list of publications to gene_list.txt file
with open("paper_list.txt","w") as file:
	for paper in paper_list:
		file.write("{}\n".format(paper.strip("'").strip('"')))
file.close()




