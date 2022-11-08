# monarch-gene-mapping

Code for mapping source namespaces to preferred namespacing

## Installation

```bash
poetry install
```

## Usage

```bash
python -m monarch_gene_mapping.main --help
```

is a simple UI for processing the mapping data. 

## Special Data Considerations

The UniProtKB ID mappings file is huge: about an eleven (11) gigabyte _gzip_ compressed archive (as of November 2022). 
The Monarch Initiative only targets a subset of species in this file. The standard procedure is to 'pre-filter' the
data after downloading but before continued processing. This is the default 'generate' process, but the
monarch_gene_mapping.main script allows for discrete processing of this step.

