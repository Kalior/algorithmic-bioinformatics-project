import numpy as np
from typing import Iterator, Tuple


def window_idx(
    length: int, window_size: int, overlap: int
) -> Iterator[Tuple[int, int]]:
    for start in np.arange(0, length, window_size - overlap):
        end = min(length, start + window_size)
        yield start, end
