#!/usr/bin/env python
from argparse import Namespace, ArgumentParser
from typing import Generator, Dict, Set
import pandas as pd

import logging
logger = logging.getLogger(__name__)

# TODO: this is hard coded; could be parameterized for
#       user access, or set to zero to ensure full capture
MAX_ENTRIES = 0
PINC: int = 10000


def get_parameters() -> Namespace:
    parser = ArgumentParser(
        prog='Analyse Monarch Dangling Edges',
        description='Analyse gene_mappings missing from a Monarch ingest',
        epilog='Copyright 2023 - the Monarch Initiative'
    )
    parser.add_argument(
        "-f", "--filename",
        help="Monarch dangling edges GZ archive file (Default: 'monarch-kg-dangling-edges.tsv.gz')"
    )
    parser.add_argument(
        '-s', '--knowledge_source',
        help="Simple partial match text filter on the 'primary knowledge source' field of the data file."
    )
    parser.add_argument(
        '-g', '--gene_filter',
        help="Simple partial match text filter on gene names (e.g. generally, the namespaces, i.e. 'ENSEMBL')."
    )

    args = parser.parse_args()
    print(
        f"File Name:\t'{args.filename}'\n",
        f"Primary Knowledge Source:\t{args.knowledge_source}\n",
        f"Gene Filter:\t{args.gene_filter}"
    )
    return args


def read_entries_by_knowledge_source(df: pd.DataFrame, knowledge_source: str) -> Generator:
    """
    Read records from pandas.DataFrame and yield records.
    :param df: pandas.DataFrame, Dataframe containing records that represent dangling edge entries.
    :param knowledge_source: str, knowledge source filter for records
    :return: A generator for dangling edge entries from specified knowledge source or all (if former is unset)
    """
    for entry in df.to_dict("records"):
        if entry and (
                not knowledge_source or
                'primary_knowledge_source' not in entry or
                (
                        entry['primary_knowledge_source'] and
                        entry['primary_knowledge_source'].find(knowledge_source) >= 0
                )
        ):
            yield entry


def parse(filename: str, knowledge_source: str) -> Generator:
    """
    This method reads from a TSV/CSV and yields records.

    :param filename: The GZ filename to parse
    :param knowledge_source: str, knowledge source filter for records
    :return: A generator for dangling edge records
    """
    file_iter = pd.read_table(
        filename,
        dtype=str,
        chunksize=10000,
        keep_default_na=False,
        compression='gzip'
    )
    for chunk in file_iter:
        yield from read_entries_by_knowledge_source(chunk, knowledge_source)


def process_entry(entry: Dict, gene_filter: str, result_set: Set[str]):
    """
    Dump selected fields of the dangling edge.

    :param entry: Dict entry to be processed
    :param gene_filter: str, filter on gene names in entry
    :param result_set: Set[str] filter on gene names in entry
    :return: None
    """
    subject_id: str = entry['subject']
    if subject_id and (not gene_filter or subject_id.find(gene_filter) >= 0):
        result_set.add(subject_id)
    object_id: str = entry['object']
    if object_id and (not gene_filter or object_id.find(gene_filter) >= 0):
        result_set.add(object_id)


def analyse_file(filename: str, knowledge_source: str, gene_filter: str) -> Set[str]:
    """
    Opens and reads a "dangling edges" TSV file.

    :param filename: str, file path to GZ archive of input data
    :param knowledge_source: str, simple string filter on primary knowledge source whose records are being extracted
    :param gene_filter: simple string filter on gene names to be compiled
    :return: Set[str] of target identifiers
    """
    n: int = 0
    entry: Dict
    result_set: Set[str] = set()
    for entry in parse(filename, knowledge_source):
        process_entry(entry, gene_filter, result_set)
        n += 1
        if not (n % PINC):
            print('.', end="" if n % (PINC*40) else "\n", flush=True)
        if MAX_ENTRIES and n >= MAX_ENTRIES:
            break

    print(f"\n{n} '{knowledge_source}' entries assessed from input '{filename}'")
    return result_set


def dump_results(knowledge_source: str, result_set: Set):
    """
    Dumps the resulting set of extracted identifiers into a text file.

    :param knowledge_source: str, primary knowledge source origin of the data
    :param result_set: Set of unmapped identifiers matching to the 'gene_filter' string
    :return:
    """
    with open(f"{knowledge_source}_unmapped_gene_ids.txt", mode="w") as output:
        for identifier in result_set:
            print(identifier, file=output)


def main():

    args = get_parameters()
    result_set: Set[str] = analyse_file(args.filename, args.knowledge_source, args.gene_filter)
    dump_results(args.knowledge_source, result_set)


if __name__ == "__main__":
    main()
