import gzip, json, shutil
import pandas as pd
from pandas.core.frame import DataFrame

def alliance_ncbi_mapping() -> DataFrame:
    alliance_files = [
         'BGI_FB.json.gz',
         'BGI_MGI.json.gz',
         'BGI_RGD.json.gz',
         'BGI_SGD.json.gz',
         'BGI_WB.json.gz',
         'BGI_XB.json.gz',
         'BGI_ZFIN.json.gz',
    ]

    data = []
    for file in alliance_files:
        with gzip.open(f"data/alliance/{file}", 'rt', encoding='UTF-8') as zipfile:
            with open(f"data/alliance/{file[:-3]}", 'wt') as f_out:
                shutil.copyfileobj(zipfile, f_out)
            # genes = json.load(zipfile)
            # data.append(genes)

    # gene_maps = pd.DataFrame(data)
    # return gene_maps

gene_maps = alliance_ncbi_mapping()
print(gene_maps)