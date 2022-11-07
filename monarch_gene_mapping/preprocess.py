"""
Utility methods for pre-filtering of huge input files of
mapping data (i.e. data/uniprot/idmapping_selected.tab.gz)
"""
from sys import stderr
from os import sep, remove
from tarfile import (
    open as tar_open,
    ReadError,
    CompressionError,
    TarInfo

)
from datetime import datetime

from io import TextIOWrapper
import logging
from typing import Optional, List

logger = logging.getLogger()


_ncbitaxon_catalog = {
    # TODO: may need to further build up this catalog
    #       of species - just starting with the STRING DB list + Dictyostelium
    "HUMAN": "9606",
    "MOUSE": "10090",
    "CANLF": "9615",    # Canis lupus familiaris - domestic dog
    # "FELCA": "9685",  # Felis catus - domestic cat
    "BOVIN": "9913",    # Bos taurus - cow
    "PIG": "9823",      # Sus scrofa - pig
    "RAT": "10116",
    "CHICK": "9031",
    "XENTR": "8364",   # Xenopus tropicalis - tropical clawed frog
    "DANRE": "7955",
    "DROME": "7227",
    "CAEEL": "6239",
    "DICDI": "44689",
    "EMENI": "227321",  # Emericella nidulans (strain FGSC A4 etc.) (Aspergillus nidulans)
    "SCHPO": "4896",
    "YEAST": "4932",
}


def target_taxon(line: Optional[str]) -> bool:
    if not line:
        return False
    part: List[str] = line.split("|")
    if part[0] in _ncbitaxon_catalog.keys():
        return True
    else:
        return False


def filter_uniprot_id_mapping_file(
        directory: str = '.',
        source_filename: str = "data/uniprot/idmapping_selected.tab.gz",
        target_filename: str = "data/uniprot/idmapping_filtered.tab.gz",
        number_of_lines: int = 0
) -> bool:
    """
    Filters contents of a Panther orthologs tar.gz archive against the target list of taxa.
    :param directory: str, location of source data file
    :param source_filename: str, source data file name
    :param target_filename: str, target data file name
    :param number_of_lines: int, number of lines parsed; 'all' lines parsed if omitted or set to zero
    :return: bool, True if filtering was successful; False if unsuccessful
    """
    assert source_filename
    assert target_filename
    assert number_of_lines >= 0
    print(
        f"\nBegin file filtering '{number_of_lines if number_of_lines else 'all'}'" +
        f" lines in '{source_filename}' at {datetime.now().isoformat()}. " +
        f"\nPatience! This may take a little awhile!...", file=stderr
    )
    # A standard TAR file with a single entry identical in name to filename
    source_tar_filename: str = f"{source_filename}.tar.gz"
    source_tar_file_path: str = f"{directory}{sep}{source_tar_filename}"
    target_file_path: str = f"{directory}{sep}{target_filename}"
    n: int = 0
    try:
        with tar_open(source_tar_file_path, mode='r') as source_tar_file:
            with open(target_file_path, mode='w', encoding='utf-8') as target_file:
                # first member assumed to be the one and only target member
                member = source_tar_file.next()
                if not member:
                    logger.error(f"filter_file() tar archive '{source_tar_filename}' is empty?")
                    print(
                        f"Failed preprocessing '{source_filename}' at {datetime.now().isoformat()}", file=stderr
                    )
                    return False
                with source_tar_file.extractfile(member) as orthologs_reader:  # io.BufferedReader
                    orthologs_file = TextIOWrapper(orthologs_reader)
                    for line in orthologs_file:
                        if target_taxon(line):
                            print(line, file=target_file, end="")
                            n += 1
                        if number_of_lines and n > number_of_lines:
                            break
    except (ReadError, CompressionError) as e:
        logger.error(f"filter_file() tar file access exception: {str(e)}")
        print(f"Failed preprocessing '{source_filename}' at {datetime.now().isoformat()}", file=stderr)
        return False

    target_tar_filename: str = f"{target_filename}.tar.gz"
    target_tar_file_path: str = f"{directory}{sep}{target_tar_filename}"
    with tar_open(target_tar_file_path, mode='w:gz') as target_tar_file:
        target_tarinfo: TarInfo = target_tar_file.gettarinfo(name=target_file_path, arcname=target_filename)
        with open(target_file_path, mode='rb') as target_file:
            target_tar_file.addfile(target_tarinfo, fileobj=target_file)

    # delete the naked file
    remove(target_file_path)

    print(f"Finished filtering file '{source_filename}' at {datetime.now().isoformat()}", file=stderr)
    return True
