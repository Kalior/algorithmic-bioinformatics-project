from typing import List

from Bio.SeqRecord import SeqRecord
from ete3 import NCBITaxa
import pandas as pd


def get_metadata(records: List[SeqRecord]):
    ncbi = NCBITaxa()

    species = [gb.annotations["organism"] for gb in records]
    name_translator = ncbi.get_name_translator(species)

    sought_ranks = ["superkingdom", "order", "family", "subfamily", "genus", "species"]

    metadata = []

    for gb in records:
        taxid = name_translator[gb.annotations["organism"]][0]
        lineage = ncbi.get_lineage(taxid)
        ranks = ncbi.get_rank(lineage)
        names = ncbi.get_taxid_translator(lineage)
        taxonomy = {ranks[k]: names[k] for k in lineage if ranks[k] in sought_ranks}
        metadata.append({**taxonomy, "aid": gb.id})

    df = pd.DataFrame(metadata)
    df.to_csv("metadata.csv")

    return df
