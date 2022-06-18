from dataclasses import dataclass

import numpy as np

from generator.cn_sampler import CNSampler


@dataclass
class CountsGenerator:
    cn_sampler: CNSampler

    def add_noise_to_counts(self, counts: np.ndarray) -> np.ndarray:
        return np.vectorize(lambda x: self.cn_sampler.sample_corrected_count(x))(counts)


