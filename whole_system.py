from Bio import Entrez

Entrez.email = "dan_egan@live.co.uk"
handle = Entrez.efetch(db="nucleotide", id="AY851612", rettype="gb", retmode="text")
print(handle.readline())
handle.close()

