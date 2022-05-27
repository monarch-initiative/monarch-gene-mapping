import gzip
import json
from typing import Dict
import typer
import pathlib

import pandas as pd
from pandas.core.frame import DataFrame

typer_app = typer.Typer()


def bgi2sssom(bgi) -> Dict:
    """
    Convert a single Alliance BGI document to a SSSOM row
    :param bgi:
    :return:
    """
    gene = bgi['basicGeneticEntity']
    xref = [x['id'].replace("NCBI_Gene", "NCBIGene").replace("NCBI_GENE", "NCBIGene")
            for x in gene['crossReferences']
            if "NCBI" in x['id']]
    prim_id = gene['primaryId'].replace("DRSC:XB:", "Xenbase:")
    if len(xref):
        ncbi_xref = xref[0]
    else:
        return None
    return {
        "subject_id": ncbi_xref,
        "predicate_id": "skos:exactMatch",
        "object_id": prim_id
    }


def alliance_ncbi_mapping() -> DataFrame:
    alliance_files = [
         'BGI_MGI.json.gz',
         'BGI_RGD.json.gz',
         'BGI_ZFIN.json.gz',
         'BGI_FB.json.gz',
         'BGI_SGD.json.gz',
         'BGI_WB.json.gz',
         # 'BGI_XB.json.gz'
    ]

    data = []
    for file in alliance_files:
        with gzip.open(f"data/alliance/{file}", 'rt', encoding='UTF-8') as zipfile:
            genes = json.load(zipfile)
        for index, bgi in enumerate(genes['data']):
            sssom = bgi2sssom(bgi)
            if sssom is not None:
                data.append(sssom)

    gene_maps = pd.DataFrame(data)
    return gene_maps


def hgnc_ncbi_mapping() -> DataFrame:
    hgnc = pd.read_csv("data/hgnc/hgnc_complete_set.txt", sep="\t", dtype="string")
    hgnc.rename(columns={"entrez_id": "subject_id", "hgnc_id": "object_id"}, inplace=True)
    hgnc["subject_id"] = 'NCBIGene:' + hgnc["subject_id"]
    hgnc["predicate"] = "skos:exactMatch"
    hgnc = hgnc[["subject_id", "predicate", "object_id"]]
    hgnc.dropna(inplace=True)

    return hgnc


def generate_gene_mappings() -> DataFrame:

    mapping_dataframes = []

    alliance_ncbi = alliance_ncbi_mapping()
    assert(len(alliance_ncbi) > 200000)
    mapping_dataframes.append(alliance_ncbi)

    hgnc_ncbi = hgnc_ncbi_mapping()
    assert(len(hgnc_ncbi) > 40000)
    mapping_dataframes.append(hgnc_ncbi)

    mappings = pd.concat(mapping_dataframes)

    return mappings


@typer_app.command()
def main(output_dir="output"):
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    mappings = generate_gene_mappings()
    mappings.to_csv(f"{output_dir}/gene_mappings.tsv", sep="\t", index=False)


if __name__ == "__main__":
    typer_app()
