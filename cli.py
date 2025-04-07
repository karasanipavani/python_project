import argparse
from pubmed_paper_fetcher.fetcher import fetch_pubmed_ids, fetch_pubmed_details
from pubmed_paper_fetcher.parser import parse_article
from pubmed_paper_fetcher.utils import write_to_csv, print_to_console


def main():
    parser = argparse.ArgumentParser(description="Fetch PubMed papers with non-academic authors.")
    parser.add_argument("query", help="PubMed query string")
    parser.add_argument("-d", "--debug", action="store_true", help="Enable debug output")
    parser.add_argument("-f", "--file", help="Filename to save CSV results")
    args = parser.parse_args()

    ids = fetch_pubmed_ids(args.query)
    if args.debug:
        print(f"Fetched {len(ids)} IDs")

    articles = fetch_pubmed_details(ids)
    results = []

    for art in articles:
        parsed = parse_article(art)
        if parsed:
            results.append(parsed)

    if args.file:
        write_to_csv(results, args.file)
    else:
        print_to_console(results)
