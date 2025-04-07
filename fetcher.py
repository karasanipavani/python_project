from typing import List
from Bio import Entrez

Entrez.email = "your.email@example.com"  # Replace with your email


def fetch_pubmed_ids(query: str, max_results: int = 50) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record["IdList"]


def fetch_pubmed_details(id_list: List[str]) -> List[dict]:
    handle = Entrez.efetch(db="pubmed", id=",".join(id_list), rettype="medline", retmode="xml")
    records = Entrez.read(handle)
    return records["PubmedArticle"]
