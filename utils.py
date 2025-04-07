from typing import List
import csv

def write_to_csv(data: List[dict], filename: str) -> None:
    if not data:
        print("No results to write.")
        return

    with open(filename, mode="w", newline="", encoding="utf-8") as file:
        writer = csv.DictWriter(file, fieldnames=data[0].keys())
        writer.writeheader()
        for row in data:
            writer.writerow(row)


def print_to_console(data: List[dict]) -> None:
    for row in data:
        print("\n---\n")
        for k, v in row.items():
            print(f"{k}: {v}")
