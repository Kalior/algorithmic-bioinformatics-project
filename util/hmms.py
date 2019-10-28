import logging
import os
import pickle
from datetime import datetime
from typing import Dict, List, Union

import numpy as np
import pandas as pd
from Bio.SeqRecord import SeqRecord

from hmmlearn import hmm

from . import rolling


def test(fasta_id: str, records: List[SeqRecord]):
    metadata = pd.read_csv("metadata.csv", index_col=0)

    now = datetime.now().strftime("%Y-%m-%d")
    output_directory = os.path.join("results", now)

    if not os.path.isdir(output_directory):
        os.makedirs(output_directory)

    model = load_hmm(fasta_id)

    print(
        {
            k: list(v.values())[0]
            for k, v in metadata[metadata.aid == fasta_id].to_dict().items()
        }
    )

    scores = []
    for record in records:
        scores.extend(rolling_score(record, model, metadata))

    df = pd.DataFrame(scores)

    output_file = os.path.join(
        output_directory, "output-{}.csv".format(datetime.now().strftime("%H-%M"))
    )
    df.to_csv(output_file)

    print("Saved to {}".format(output_file))

    return df


def rolling_score(
    record: SeqRecord,
    model: hmm.MultinomialHMM,
    metadata: pd.DataFrame,
    window_size: int = 1000,
    overlap: int = 950,
) -> List[Dict[str, Union[str, float, int]]]:
    scores = []

    enc = {"A": 0, "C": 1, "G": 2, "T": 3}
    sequence = np.array([enc.get(c, 0) for c in str(record.seq)])

    for start, end in rolling.window_idx(len(sequence), window_size, overlap):
        subsequence = sequence[start:end].reshape(-1, 1)

        score = model.score(subsequence) / (end - start)
        scores.append(
            {
                "id": record.id,
                "start": start,
                "score": score,
                "relative_start": start / len(sequence),
                **{
                    k: list(v.values())[0]
                    for k, v in metadata[metadata.aid == record.id].to_dict().items()
                },
            }
        )

    return scores


def train_hmms(records: List[SeqRecord]) -> List[hmm.MultinomialHMM]:
    enc = {"A": 0, "C": 1, "G": 2, "T": 3}

    models = []
    hmms_directory = "../hmms"
    for record in records:
        logging.info("Training HMM for {}".format(record.id))
        sequence = np.array([enc.get(c, 0) for c in str(record.seq)]).reshape(-1, 1)

        model = hmm.MultinomialHMM(n_components=10, n_iter=100).fit(sequence)
        models.append(model)
        with open(os.path.join(hmms_directory, record.id + ".pkl"), "wb") as file_:
            pickle.dump(model, file_)

    return models


def load_hmm(fasta_id: str, hmms_directory: str = "../hmms") -> hmm.MultinomialHMM:
    with open(os.path.join(hmms_directory, fasta_id + ".pkl"), "rb") as file_:
        model = pickle.load(file_)

    return model
