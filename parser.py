from typing import List, Tuple, Optional
import re

NON_ACADEMIC_KEYWORDS = ["Inc", "Ltd", "LLC", "Pharma", "Biotech", "Corporation", "Corp"]
ACADEMIC_KEYWORDS = ["University", "College", "Institute", "Hospital", "School", "Center", "Centre"]


def is_non_academic(affiliation: str) -> bool:
    return (
        any(kw in affiliation for kw in NON_ACADEMIC_KEYWORDS)
        and not any(kw in affiliation for kw in ACADEMIC_KEYWORDS)
    )


def parse_article(article: dict) -> Optional[dict]:
    try:
        medline = article["MedlineCitation"]
        article_data = medline["Article"]
        pmid = str(medline["PMID"])
        title = article_data.get("ArticleTitle", "")
        pub_date = article_data.get("Journal", {}).get("JournalIssue", {}).get("PubDate", {})
        pub_date_str = pub_date.get("Year", "")

        authors = article_data.get("AuthorList", [])
        non_academic_authors = []
        company_affiliations = []
        corresponding_email = None

        for author in authors:
            if "AffiliationInfo" in author:
                for aff in author["AffiliationInfo"]:
                    affiliation = aff.get("Affiliation", "")
                    if is_non_academic(affiliation):
                        name = author.get("LastName", "") + ", " + author.get("ForeName", "")
                        non_academic_authors.append(name)
                        company_affiliations.append(affiliation)
                        email_match = re.search(r"[\w\.-]+@[\w\.-]+", affiliation)
                        if email_match:
                            corresponding_email = email_match.group(0)

        if not non_academic_authors:
            return None

        return {
            "PubmedID": pmid,
            "Title": title,
            "Publication Date": pub_date_str,
            "Non-academic Author(s)": "; ".join(non_academic_authors),
            "Company Affiliation(s)": "; ".join(set(company_affiliations)),
            "Corresponding Author Email": corresponding_email or "N/A",
        }
    except Exception:
        return None
