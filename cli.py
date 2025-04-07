# pubmed_paper_fetcher/cli.py

import argparse
import csv
from pubmed_paper_fetcher.fetcher import fetch_pubmed_ids, fetch_pubmed_details, parse_author_data
from typing import List
import logging

def write_csv(data: List[dict], filename: str):
    keys = ["PubmedID", "Title", "Publication Date", "Non-academic Author(s)", "Company Affiliation(s)", "Corresponding Author Email"]
    with open(filename, "w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=keys)
        writer.writeheader()
        for item in data:
            item["Non-academic Author(s)"] = ", ".join(item["Non-academic Author(s)"])
            item["Company Affiliation(s)"] = ", ".join(item["Company Affiliation(s)"])
            writer.writerow(item)

def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers with non-academic authors from pharma/biotech companies.")
    parser.add_argument("query", type=str, help="PubMed search query")
    parser.add_argument("-f", "--file", type=str, help="Output CSV filename")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug logging")

    args = parser.parse_args()

    if args.debug:
        logging.basicConfig(level=logging.DEBUG)

    ids = fetch_pubmed_ids(args.query)
    logging.info(f"Found {len(ids)} results")

    results = []
    for pmid in ids:
        raw_record = fetch_pubmed_details(pmid)
        if raw_record:
            parsed = parse_author_data(raw_record)
            if parsed["Company Affiliation(s)"]:
                results.append(parsed)

    if args.file:
        write_csv(results, args.file)
        print(f"Saved {len(results)} results to {args.file}")
    else:
        for r in results:
            print(r)
