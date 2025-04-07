# pubmed_paper_fetcher/fetcher.py

from typing import List, Dict, Optional
from Bio import Entrez
import re
import logging

Entrez.email = "pavanikar@gmail.com"

NON_ACADEMIC_KEYWORDS = ["pharma", "biotech", "inc", "ltd", "gmbh", "corp", "llc"]

def is_non_academic(affiliation: str) -> bool:
    university_keywords = ["university", "college", "institute", "school", "hospital", "center", "centre", "lab"]
    return not any(keyword.lower() in affiliation.lower() for keyword in university_keywords)

def is_pharma_or_biotech(affiliation: str) -> bool:
    return any(keyword in affiliation.lower() for keyword in NON_ACADEMIC_KEYWORDS)

def fetch_pubmed_ids(query: str, max_results: int = 100) -> List[str]:
    handle = Entrez.esearch(db="pubmed", term=query, retmax=max_results)
    record = Entrez.read(handle)
    return record.get("IdList", [])

def fetch_pubmed_details(pmid: str) -> Optional[Dict]:
    handle = Entrez.efetch(db="pubmed", id=pmid, rettype="medline", retmode="text")
    record = handle.read()
    handle.close()
    return record

def parse_author_data(record: str) -> Dict:
    from Bio import Medline
    from io import StringIO

    data = list(Medline.parse(StringIO(record)))[0]

    result = {
        "PubmedID": data.get("PMID", ""),
        "Title": data.get("TI", ""),
        "Publication Date": data.get("DP", ""),
        "Non-academic Author(s)": [],
        "Company Affiliation(s)": [],
        "Corresponding Author Email": ""
    }

    authors = data.get("AU", [])
    affiliations = data.get("AD", [])
    emails = re.findall(r'[\w\.-]+@[\w\.-]+', affiliations) if isinstance(affiliations, str) else []

    # Extract useful info
    if isinstance(affiliations, str):
        affiliations = [affiliations]

    for author, affiliation in zip(authors, affiliations):
        if is_non_academic(affiliation):
            result["Non-academic Author(s)"].append(author)
            if is_pharma_or_biotech(affiliation):
                result["Company Affiliation(s)"].append(affiliation)

    if emails:
        result["Corresponding Author Email"] = emails[0]

    return result
