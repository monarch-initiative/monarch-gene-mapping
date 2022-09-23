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
            if x['id'].startswith("NCBI_Gene:")
            or x['id'].startswith("UniProtKB:")
            or x['id'].startswith("ENSEMBL:")]
    prim_id = gene['primaryId'].replace("DRSC:XB:", "Xenbase:")
    if len(xref):
        ncbi_xref = xref[0]
    else:
        return None
    return {
        "subject_id": prim_id,
        "predicate_id": "skos:exactMatch",
        "object_id": ncbi_xref,
        "mapping_justification": "semapv:UnspecifiedMatching",
    }


def alliance_ncbi_mapping() -> DataFrame:
    alliance_files = [
        'BGI_MGI.json.gz',
        'BGI_RGD.json.gz',
        'BGI_ZFIN.json.gz',
        'BGI_FB.json.gz',
        'BGI_SGD.json.gz',
        'BGI_WB.json.gz',
        'BGI_XB.json.gz'
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


def hgnc_mapping(subject_column,
                 object_curie_prefix,
                 predicate_id="skos:exactMatch",
                 object_list_delimiter=None) -> DataFrame:
    hgnc = pd.read_csv("data/hgnc/hgnc_complete_set.txt", sep="\t", dtype="string")

    hgnc.rename(columns={"hgnc_id": "subject_id", subject_column: "object_id"}, inplace=True)

    # if the object column has a list, use explode to roll it down and make a SSSOM triple for each
    if object_list_delimiter is not None \
            and len(hgnc[hgnc['object_id'].str.contains(object_list_delimiter)]) > 0:
        hgnc = hgnc.assign(object_id=hgnc.object_id.str.split("|")).explode('object_id')

    hgnc["object_id"] = object_curie_prefix + hgnc["object_id"]
    hgnc["predicate_id"] = predicate_id
    hgnc["mapping_justification"] = "semapv:UnspecifiedMatching"
    hgnc = hgnc[["subject_id", "predicate_id", "object_id", "mapping_justification"]]
    hgnc.dropna(inplace=True)

    return hgnc


def generate_gene_mappings() -> DataFrame:
    mapping_dataframes = []

    ncbi_to_alliance = alliance_ncbi_mapping()
    assert (len(ncbi_to_alliance) > 200000)
    mapping_dataframes.append(ncbi_to_alliance)

    ncbi_to_hgnc = hgnc_mapping("entrez_id", "NCBIGene:")
    assert (len(ncbi_to_hgnc) > 40000)
    mapping_dataframes.append(ncbi_to_hgnc)

    omim_to_hgnc = hgnc_mapping("omim_id", "OMIM:", object_list_delimiter="\\|")
    assert (len(omim_to_hgnc) > 16000)
    mapping_dataframes.append(omim_to_hgnc)

    uniprot_to_hgnc = hgnc_mapping("uniprot_ids",
                                   "UniProtKB:",
                                   predicate_id="skos:closeMatch",
                                   object_list_delimiter="\\|")
    assert (len(uniprot_to_hgnc) > 20000)
    mapping_dataframes.append(uniprot_to_hgnc)

    mappings = pd.concat(mapping_dataframes)

    return mappings


@typer_app.command()
def main(output_dir="output"):
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    mappings = generate_gene_mappings()
    mappings.to_csv(f"{output_dir}/gene_mappings.tsv", sep="\t", index=False)


if __name__ == "__main__":
    typer_app()
