from typing import List, Dict, Union

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord

from . import rolling


def rolling_gc(
    record: SeqRecord,
    metadata: pd.DataFrame,
    window_size: int = 1000,
    overlap: int = 950,
) -> List[Dict[str, Union[str, float, int]]]:
    gc_conversion = {"A": 0, "C": 1, "G": 1, "T": 0}
    gc_sequence = np.array([gc_conversion.get(c, 0) for c in record.seq])
    gc_mean = np.mean(gc_sequence)
    scores = []
    for start, end in rolling.window_idx(len(gc_sequence), window_size, overlap):
        subsequence = gc_sequence[start:end]

        score = np.mean(subsequence)
        scores.append(
            {
                "id": record.id,
                "start": start,
                "score": score - gc_mean,
                "raw_score": score,
                **{
                    k: list(v.values())[0]
                    for k, v in metadata[metadata.aid == record.id].to_dict().items()
                },
            }
        )

    return scores
