#!/bin/env python
import sys
from argparse import Namespace, ArgumentParser
from typing import List


MAX_LINES = 0


def get_parameters() -> Namespace:
    print("Entering 'get_parameters()'")
    parser = ArgumentParser(
        prog='analyse mappings',
        description='Analyses dangling edges from a Monarch ingest',
        epilog='Copyright 2023 - the Monarch Initiative'
    )
    parser.add_argument('filename')
    parser.add_argument('-s', '--primary_knowledge_source')
    parser.add_argument('-v', '--verbose', action='store_true')  # on/off flag

    args = parser.parse_args()
    return args


def valid_entry(data: List[str]) -> bool:
    # simple check for empty or null data fields
    # if all([bool(item) for item in data]):
    # print(data[1], data[3])
    if data[1] and data[3]:
        return True
    else:
        return False


def dump_entry(fields: List[str], line: str):
    # print(line)
    data: List[str] = line.split("\t")
    if not valid_entry(data):
        for f in range(0, len(data)):
            value: str = data[f]  # .replace(':', '\\:')
            print(f"{fields[f]:>30}: {value}")
        print()
        return True
    return False


def analyse_file(filename: str, primary_knowledge_source: str):
    print(f"Entering 'analyse_file({filename} filtering on {primary_knowledge_source})'")
    n: int = 0
    d: int = 0
    # with open(sys.stdin) as input_file:
    headings: str = sys.stdin.readline()
    fields: List[str] = [tag.strip() for tag in headings.split("\t")]
    for line in sys.stdin:
        if not line.find(primary_knowledge_source):
            continue
        n += 1
        if dump_entry(fields, line):
            d += 1
        if MAX_LINES and n > MAX_LINES:
            break
    print(f"Found '{d}' invalid entries in '{n}' lines assessed from '{filename}'")


def main():
    print("Entering 'main()'")
    args = get_parameters()
    print(
        f"File Name:\t'{args.filename}'\n",
        f"Primary KS:\t{args.primary_knowledge_source}",
        f"Verbose?\t'{args.verbose}'"
    )
    analyse_file(args.filename, args.primary_knowledge_source)


if __name__ == "__main__":
    main()

