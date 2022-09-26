import typer, pathlib
from kghub_downloader.download_utils import download_from_yaml
from monarch_gene_mapping.cli_utils import * 

import logging
typer_app = typer.Typer()


@typer_app.command(name="download-data")
def _download():
    download_from_yaml(
        yaml_file='monarch_gene_mapping/download.yaml',
        output_dir='.',
    )

@typer_app.command()
def generate(
    output_dir = typer.Option("output", help="Output directory"),
    download: bool = typer.Option(False, help="Pass to first download required data ")
    ):
    if download:
        _download()
        print("\nData download complete!\n")
    print("\nGenerating gene mapping...\n")
    pathlib.Path(output_dir).mkdir(parents=True, exist_ok=True)
    mappings = generate_gene_mappings()
    mappings.to_csv(f"{output_dir}/gene_mappings.tsv", sep="\t", index=False)


if __name__ == "__main__":
    typer_app()
