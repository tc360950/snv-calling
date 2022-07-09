from dataclasses import dataclass


@dataclass
class Parameters:
    m: float  # per allele coverage
    e: float  # sequencing error
    q: float  # read success probablity
