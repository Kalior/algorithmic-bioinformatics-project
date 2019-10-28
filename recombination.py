from util import kmers, metadata

from Bio import SeqIO
import scipy
import pandas as pd
import numpy as np
import seaborn as sns
import matplotlib.pyplot as plt

gb_path = "herpersviridae.gb"
records = {gb.id: gb for gb in SeqIO.parse(gb_path, "genbank")}
metadata = metadata.get_metadata(records.values())
scores = kmers.rolling_k_mer_composition(
    records["NC_007605.1"], records["NC_014567.1"], 3, metadata
)

tetramer = pd.DataFrame(scores)
one = np.array([v for v in tetramer[tetramer.type == "NC_007605.1"]["kmers"].values])
two = np.array([v for v in tetramer[tetramer.type == "NC_014567.1"]["kmers"].values])

_, p_values = scipy.stats.ttest_ind(one, two, axis=1)
print(p_values)
p_values = sorted(p_values)
print(p_values)

df = tetramer.loc[tetramer.type == "NC_007605.1", :]
df["p_values"] = p_values
sns.lineplot(x="start", y="p_values", data=df)
