import gzip, json
from typing import Dict
import pandas as pd
from pandas.core.frame import DataFrame
from typing import List


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
        'BGI_FB.json.gz',
        'BGI_MGI.json.gz',
        'BGI_RGD.json.gz',
        'BGI_SGD.json.gz',
        'BGI_WB.json.gz',
        'BGI_XBXT.json.gz',
        'BGI_XBXL.json.gz',
        'BGI_ZFIN.json.gz',
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


def add_prefix(prefix, column) -> pd.Series:
    return prefix + column.astype('str')


def ensemble_ncbi_mappings(
        subject_column: str = "GeneID",
        subject_curie_prefix: str = "NCBIGene:",
        object_column: str = "Ensembl_gene_identifier",
        object_curie_prefix: str = "ENSEMBL:",
        predicate_id: str = "skos:exactMatch",
        mapping_justification: str = "semapv:UnspecifiedMatching",
        tax_ids: List[int] = [9031]) -> DataFrame:

    ensembl_2_gene = pd.read_csv("data/ensembl/gene2ensembl.gz", compression="gzip", sep="\t")
    ensembl_filtered = ensembl_2_gene.loc[ensembl_2_gene["#tax_id"].isin(tax_ids)].copy()

    ensembl_filtered.loc[:, subject_column] = add_prefix(subject_curie_prefix, ensembl_filtered[subject_column])
    ensembl_filtered.loc[:, object_column] = add_prefix(object_curie_prefix, ensembl_filtered[object_column])
    ensembl_filtered.loc[:, "predicate_id"] = predicate_id
    ensembl_filtered.loc[:, "mapping_justification"] = mapping_justification

    columns = {subject_column: "subject_id", object_column: "object_id"}
    select_columns = ["subject_id", "predicate_id", "object_id", "mapping_justification"]
    ensembl_mappings = ensembl_filtered.rename(columns=columns).loc[:, select_columns].drop_duplicates()

    return ensembl_mappings


def generate_gene_mappings() -> DataFrame:

    mapping_dataframes = []

    ncbi_to_alliance = alliance_ncbi_mapping()
    assert(len(ncbi_to_alliance) > 200000)
    mapping_dataframes.append(ncbi_to_alliance)

    ncbi_to_hgnc = hgnc_mapping("entrez_id", "NCBIGene:")
    assert(len(ncbi_to_hgnc) > 40000)
    mapping_dataframes.append(ncbi_to_hgnc)

    omim_to_hgnc = hgnc_mapping("omim_id", "OMIM:", object_list_delimiter="\\|")
    assert(len(omim_to_hgnc) > 16000)
    mapping_dataframes.append(omim_to_hgnc)

    uniprot_to_hgnc = hgnc_mapping("uniprot_ids",
                                   "UniProtKB:",
                                   predicate_id="skos:closeMatch",
                                   object_list_delimiter="\\|")
    assert(len(uniprot_to_hgnc) > 20000)
    mapping_dataframes.append(uniprot_to_hgnc)

    ncbi_to_ensembl = ensemble_ncbi_mappings(
        subject_column="Ensembl_gene_identifier",
        subject_curie_prefix="ENSEMBL:",
        object_column="GeneID",
        object_curie_prefix="NCBIGene:",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching",
        tax_ids=[9031]
    )
    assert (len(ncbi_to_ensembl) > 15000)
    mapping_dataframes.append(ncbi_to_ensembl)

    mappings = pd.concat(mapping_dataframes)

    return mappings
