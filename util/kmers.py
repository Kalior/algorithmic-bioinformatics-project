import itertools
from collections import defaultdict
from functools import lru_cache

import networkx as nx
import numpy as np

from . import rolling

from Bio.SeqRecord import SeqRecord
import pandas as pd

from typing import List, Dict, Union, Tuple


def rolling_k_mer_composition(
    record: SeqRecord,
    to: SeqRecord,
    k: int,
    metadata: pd.DataFrame,
    window_size: int = 1000,
    overlap: int = 950,
) -> List[Dict[str, Union[str, float, int]]]:

    to_kmer_frequency = kmer_frequency(to.seq, k)
    all_scores = []
    for start, end in rolling.window_idx(len(record.seq), window_size, overlap):
        record_kmer_frequency = kmer_frequency(record.seq[start:end], k)

        scores = {
            "diff": to_kmer_frequency - record_kmer_frequency,
            record.id: record_kmer_frequency,
            to.id: to_kmer_frequency,
        }

        all_scores.extend(
            [
                {
                    "id": record.id,
                    "start": start,
                    "score": np.linalg.norm(kmer_freq),
                    "kmers": kmer_freq,
                    "type": type_,
                    **{
                        k: list(v.values())[0]
                        for k, v in metadata[metadata.aid == record.id]
                        .to_dict()
                        .items()
                    },
                }
                for type_, kmer_freq in scores.items()
            ]
        )

    return all_scores


def kmer_frequency(sequence: str, k: int, normalise: bool = True) -> np.ndarray:
    mono_counts = defaultdict(int)
    kmer_counts = defaultdict(int)

    alphabet = ["A", "C", "G", "T"]

    sequence_length = len(sequence)

    for i in range(sequence_length - k):
        end_i = i + k
        kmer = sequence[i:end_i]
        if all(c in alphabet for c in kmer):
            kmer_counts[kmer] += 1

    for c in sequence:
        if c in alphabet:
            mono_counts[c] += 1

    mono_frequencies = defaultdict(
        float, {key: value / sequence_length for key, value in mono_counts.items()}
    )

    kmers = ["".join(s) for s in itertools.product("ACGT", repeat=k)]
    for kmer in kmers:
        kmer_counts[kmer] += 0

    if normalise:
        return np.array(
            [
                value
                / (sequence_length - k)
                / np.prod([mono_frequencies[c] for c in key])
                for key, value in kmer_counts.items()
            ]
        )
    else:
        return np.array(
            [value / (sequence_length - k) for key, value in kmer_counts.items()]
        )


@lru_cache(None)
def get_kmers(sequence: str, k: int) -> List[Tuple[str, int]]:
    kmer_counts = defaultdict(int)

    for i in range(len(sequence) - k):
        end_i = i + k
        kmer = sequence[i:end_i]
        kmer_counts[kmer] += 1

    return list(kmer_counts.items())


def shared_kmer_network(records: List[SeqRecord], k: int) -> nx.Graph:
    graph = nx.Graph()
    graph.add_nodes_from([record.id for record in records])

    for left, right in itertools.product(records, repeat=2):
        if left.id == right.id:
            continue

        left_kmers = get_kmers(str(left.seq), k)
        right_kmers = get_kmers(str(right.seq), k)
        shared_kmers = set(kmer for kmer, _ in left_kmers).intersection(
            set(kmer for kmer, _ in right_kmers)
        )

        for kmer in shared_kmers:
            graph.add_edge(left.id, right.id, kmer=kmer)

    return graph
