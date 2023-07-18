#!/bin/env python
from argparse import Namespace, ArgumentParser
from typing import List, Generator, Optional,  Dict, Tuple
import pandas as pd

import logging
logger = logging.getLogger(__name__)

MAX_ENTRIES = 10


def get_parameters() -> Namespace:
    print("Entering 'get_parameters()'")
    parser = ArgumentParser(
        prog='Analyse Monarch Dangling Edges',
        description='Analyse gene_mappings missing from a Monarch ingest',
        epilog='Copyright 2023 - the Monarch Initiative'
    )
    parser.add_argument(
        "filename",
        default="monarch-kg-dangling-edges.tsv.gz",
        help="Monarch angling edges GZ archive file (Default: 'monarch-kg-dangling-edges.tsv.gz')"
    )
    parser.add_argument('-s', '--knowledge_source')

    args = parser.parse_args()
    return args


report_fields: Optional[List[str]] = ["subject", "object", "has_evidence"]


def read_entry(entry: Dict, knowledge_source: str) -> Optional[Dict]:
    """
    Prepare a dangling edge line record.

    :param entry: Dict, a single dangling edge record.
    :param knowledge_source: str, knowledge source filter for records
    :return: Optional[Tuple[str, Dict]], a tuple that contains node id and node data
    """
    if entry:
        # if not obj.find(knowledge_source):
        #     continue
        # node = self.validate_node(node)
        # if node:
        #     # if not None, assumed to have an "id" here...
        #     node_data = sanitize_import(node.copy(), self.list_delimiter)
        #     n = node_data["id"]
        #
        #     self.set_node_provenance(node_data)  # this method adds provided_by to the node properties/node data
        #     self.node_properties.update(list(node_data.keys()))
        #     if self.check_node_filter(node_data):
        #         self.node_properties.update(node_data.keys())
        #         return n, node_data
        if entry['primary_knowledge_source'] and entry['primary_knowledge_source'].find(knowledge_source):
            return entry
    return None


def read_entries(df: pd.DataFrame, knowledge_source: str) -> Generator:
    """
    Read records from pandas.DataFrame and yield records.
    :param df: pandas.DataFrame, Dataframe containing records that represent dangling edge entries.
    :param knowledge_source: str, knowledge source filter for records
    :return: A generator for dangling edge entries
    """
    for obj in df.to_dict("records"):
        entry: Optional[Dict] = read_entry(obj, knowledge_source)
        if entry is None:
            continue
        yield entry


def parse(filename: str, knowledge_source: str) -> Generator:
    """
    This method reads from a TSV/CSV and yields records.

    :param filename: The GZ filename to parse
    :param knowledge_source: str, knowledge source filter for records
    :return: A generator for dangling edge records
    """
    global report_fields
    file_iter = pd.read_table(
        filename,
        dtype=str,
        chunksize=10000,
        keep_default_na=False,
        compression='gzip'
    )
    for chunk in file_iter:
        # TODO: how can I capture the column identities in an ordered fashion?
        if report_fields is None:
            report_fields = list(chunk.columns)
        yield from read_entries(chunk, knowledge_source)


def dump_entry(entry: Optional[Dict], output):
    """
    Dump selected fields of the dangling edge.

    :param entry: entry to be reported
    :param output: text file device open for output
    :return:
    """
    global report_fields
    print('\t'.join([entry[f] for f in report_fields]), file=output)


def analyse_file(filename: str, knowledge_source: str, output):
    """
    Opens and reads a "dangling edges" TSV file.

    :param filename: (gzip'd) input file name
    :param knowledge_source: primary knowledge source whose records are being extracted
    :param output: text file device open for output
    :return: 
    """
    global report_fields
    print(f"Entering analyse_file({filename} filtering on {knowledge_source})")
    
    n: int = 0

    entry: Dict
    print('\t'.join(report_fields), file=output)
    for entry in parse(filename, knowledge_source):
        dump_entry(entry, output)
        n += 1
        if MAX_ENTRIES and n > MAX_ENTRIES:
            break

    print(f"Found '{n}' lines assessed from '{filename}'")


def main():

    print("Entering 'main()'")
    args = get_parameters()
    print(
        f"File Name:\t'{args.filename}'\n",
        f"Primary KS:\t{args.knowledge_source}"
    )

    with open(f"{args.knowledge_source}_dangling_edges.tsv", mode="w") as output:
        analyse_file(args.filename, args.knowledge_source, output)


if __name__ == "__main__":
    main()
