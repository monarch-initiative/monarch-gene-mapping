from os import sep
import typer
import pathlib
from kghub_downloader.download_utils import download_from_yaml
from monarch_gene_mapping.cli_utils import *
from .uniprot_idmapping_preprocess import filter_uniprot_id_mapping_file

typer_app = typer.Typer()


@typer_app.command(name="download-data")
def _download():
    download_from_yaml(
        yaml_file='monarch_gene_mapping/download.yaml',
        output_dir='.',
    )


@typer_app.command(name="filter-uniprot")
def preprocess_uniprot_mappings(
        directory: str = typer.Option(f"..{sep}data{sep}uniprot", help="Data Directory"),
        source_filename: str = typer.Option("idmapping_selected.tab", help="Target File Name"),
        target_filename: str = typer.Option("idmapping.tsv", help="Target File Name"),
        number_of_lines: int = typer.Option(0, help="Number of Lines")
):
    filter_uniprot_id_mapping_file(
        directory=directory,
        source_filename=source_filename,
        target_filename=target_filename,
        number_of_lines=number_of_lines
    )


@typer_app.command()
def generate(
        output_dir=typer.Option("output", help="Output directory"),
        download: bool = typer.Option(False, help="Pass to first download required data"),
        filter_uniprot: bool = typer.Option(True, help="Filter out UniProt ID mapping data after download")
):
    if download:
        _download()
        print("\nData download complete!\n")

    # prefilter 'target' taxa in Uniprot data
    if filter_uniprot:
        filter_uniprot_id_mapping_file(
            directory="data/uniprot",
            source_filename="idmapping_selected.tab",
            target_filename="idmapping.tsv",
            number_of_lines=0
        )

    print("\nGenerating gene mapping...\n")
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    mappings = generate_gene_mappings()

    print(f"Results saved in {output_dir}/gene_mappings.tsv")
    mappings.to_csv(f"{output_dir}/gene_mappings.tsv", sep="\t", index=False)


if __name__ == "__main__":
    typer_app()
