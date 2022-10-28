import gzip, json
import pandas as pd
from pandas.core.frame import DataFrame
from typing import Dict, List


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


def add_prefix(prefix: str, column: pd.Series) -> pd.Series:
    """
    Add a prefix to all values in a series
    :param prefix: Prefix for values in series
    :param column: Series of values for prefixing
    :return:
    """
    return prefix + column.astype('str')


def df_mappings(
        df: DataFrame,
        subject_column: str,  # "GeneID"
        object_column: str,  # "Ensembl_gene_identifier"
        predicate_id: str = "skos:exactMatch",
        mapping_justification: str = "semapv:UnspecifiedMatching",
        filter_column: str = "",  # "#tax_id",
        subject_curie_prefix: str = None,  # "NCBIGene:"
        object_curie_prefix: str = None,  # "ENSEMBL:"
        filter_ids: List[int] = None  # [9031]
) -> DataFrame:
    """
    Create specified mappings from DataFrame
    :param df: DataFrame source for mapping values
    :param subject_column: Column containing subject ID's
    :param object_column: Column containing object ID's
    :param predicate_id: String for predicate ID
    :param mapping_justification: String for mapping justification
    :param filter_column: Column to filter DataFrame on
    :param filter_ids:
    :param subject_curie_prefix: Optional curie for prefixing subject ID's
    :param object_curie_prefix: Optional curie for prefixing object ID's
    :return:
    """
    # Filtering could be extracted and done before passing df but I think there is value keeping it here.
    if filter_column != "" and isinstance(filter_ids, list):
        # Create a copy to guarantee we aren't working on a slice
        df_filtered = df.loc[df[filter_column].isin(filter_ids), :].copy()
    else:
        # Create copy so we don't modify the original DataFrame
        df_filtered = df.copy()

    df_filtered.loc[:, "predicate_id"] = predicate_id
    df_filtered.loc[:, "mapping_justification"] = mapping_justification

    columns = {subject_column: "subject_id", object_column: "object_id"}
    select_columns = ["subject_id", "predicate_id", "object_id", "mapping_justification"]
    df_select = df_filtered.rename(columns=columns).loc[:, select_columns].copy()

    if subject_curie_prefix is not None:
        df_select["subject_id"] = add_prefix(subject_curie_prefix, df_select["subject_id"])
    if object_curie_prefix is not None:
        df_select["object_id"] = add_prefix(object_curie_prefix, df_select["object_id"])
    df_map = df_select.drop_duplicates().dropna()

    return df_map


def id_list_to_sssom(df: DataFrame, column: str, delimiter: str) -> DataFrame:
    """
    Expand columns with delimiter separated lists to SSSOM triple for each
    :param df: DataFrame for column expansion
    :param column: Column name for expansion
    :param delimiter: Delimiter for splitting column
    :return:
    """
    assign_kwargs = {column: df[column].str.split(delimiter)}
    df_sssom = df.assign(**assign_kwargs).explode(column).copy()
    return df_sssom


def preprocess_alliance_df(df: DataFrame, exclude_taxon: List, include_curie: List) -> DataFrame:
    taxon_filter = ~df["TaxonID"].isin(exclude_taxon)
    curie_filter = df["GeneID"].str.contains('|'.join(include_curie))
    self_filter = df["GeneID"] != df["GlobalCrossReferenceID"]

    df_filtered = df.loc[taxon_filter & curie_filter & self_filter, :]
    return df_filtered.copy()


def generate_gene_mappings() -> DataFrame:

    mapping_dataframes = []

    ncbi_to_alliance = alliance_ncbi_mapping()
    assert(len(ncbi_to_alliance) > 200000)
    mapping_dataframes.append(ncbi_to_alliance)

    alliance_file = "data/alliance/GENECROSSREFERENCE_COMBINED.tsv.gz"
    alliance_df = pd.read_csv(alliance_file, sep="\t", dtype="string", comment='#')
    alliance_df_filtered = preprocess_alliance_df(
        df=alliance_df,
        exclude_taxon=["NCBITaxon:9606", "NCBITaxon:2697049"],
        include_curie=["MGI:", "RGD:", "FB:", "WB:", "ZFIN:", "Xenbase:"])
    ncbi_to_alliance = df_mappings(
        df=alliance_df_filtered,
        subject_column="GeneID",
        subject_curie_prefix="",
        object_column="GlobalCrossReferenceID",
        object_curie_prefix="",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching"
    )
    # assert(len(ncbi_to_alliance) > 200000)
    # mapping_dataframes.append(ncbi_to_alliance)

    hgnc_df = pd.read_csv("data/hgnc/hgnc_complete_set.txt", sep="\t", dtype="string")
    ncbi_to_hgnc = df_mappings(
        df=hgnc_df,
        subject_column="hgnc_id",
        object_column="entrez_id",
        object_curie_prefix="NCBIGene:",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching"
    )
    assert(len(ncbi_to_hgnc) > 40000)
    mapping_dataframes.append(ncbi_to_hgnc)

    omim_to_hgnc = df_mappings(
        df=id_list_to_sssom(hgnc_df, "omim_id", "|"),
        subject_column="hgnc_id",
        object_column="omim_id",
        object_curie_prefix="OMIM:",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching"
    )
    assert(len(omim_to_hgnc) > 16000)
    mapping_dataframes.append(omim_to_hgnc)

    uniprot_to_hgnc = df_mappings(
        df=id_list_to_sssom(hgnc_df, "uniprot_ids", "|"),
        subject_column="hgnc_id",
        object_column="uniprot_ids",
        object_curie_prefix="UniProtKB:",
        predicate_id="skos:closeMatch",
        mapping_justification="semapv:UnspecifiedMatching"
    )
    assert(len(uniprot_to_hgnc) > 20000)
    mapping_dataframes.append(uniprot_to_hgnc)

    ensembl_to_hgnc = df_mappings(
        df=hgnc_df,
        subject_column="hgnc_id",
        object_column="ensembl_gene_id",
        object_curie_prefix="ENSEMBL:",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching"
    )
    assert(len(ensembl_to_hgnc) > 40000)
    mapping_dataframes.append(ensembl_to_hgnc)

    ensembl_df = pd.read_csv("data/ncbi/gene2ensembl.gz", compression="gzip", sep="\t")
    ncbi_to_ensembl = df_mappings(
        df=ensembl_df,
        subject_column="GeneID",
        subject_curie_prefix="NCBIGene:",
        object_column="Ensembl_gene_identifier",
        object_curie_prefix="ENSEMBL:",
        predicate_id="skos:exactMatch",
        mapping_justification="semapv:UnspecifiedMatching",
        filter_column="#tax_id",
        filter_ids=[9031, 9615, 9913, 9823]  # Chicken: 9031, Dog: 9615, Cow, 9913, Pig: 9823
    )
    assert (len(ncbi_to_ensembl) > 70000)
    mapping_dataframes.append(ncbi_to_ensembl)


    mappings = pd.concat(mapping_dataframes)

    return mappings
