import math
import warnings
from typing import Callable

import numpy as np

from core.model_data import Parameters, EventTreeWithCounts, CellsData, create_cn_matrix


class MMEstimator:
    DEFAULT_SEQUENCING_ERROR: float = 0.0001

    @staticmethod
    def estimate(cells: CellsData, tree: EventTreeWithCounts) -> Parameters:
        CN_profiles = create_cn_matrix(cells, tree)
        zero_copy_bins = MMEstimator.__count_bin_cell_pairs_with_zero_copy(cells, CN_profiles)
        cn_sum = MMEstimator.__sum_copy_numbers(cells, CN_profiles)

        if zero_copy_bins == 0 or np.sum(cells.d[CN_profiles == 0]):
            warnings.warn(f"Can't estimate sequencing error - there are no reads for full deletions in the dataset!"
                          f"Setting sequencing error to {MMEstimator.DEFAULT_SEQUENCING_ERROR}.")
            sequencing_error = MMEstimator.DEFAULT_SEQUENCING_ERROR
        else:
            sequencing_error = np.sum(cells.d[CN_profiles == 0]) / zero_copy_bins
        per_allele_coverage = np.sum(cells.d[CN_profiles != 0]) / cn_sum
        return Parameters(
            m=per_allele_coverage,
            e=sequencing_error,
            q=0.5
        )

    @staticmethod
    def __sum_copy_numbers(cells: CellsData, CN_profiles: np.ndarray) -> int:
        return sum([cells.cell_cluster_sizes[c] * np.sum(CN_profiles[c,]) for c in range(0, cells.d.shape[0])])

    @staticmethod
    def __count_bin_cell_pairs_with_zero_copy(cells: CellsData, CN_profiles: np.ndarray) -> int:
        return sum([cells.cell_cluster_sizes[c] * np.sum(CN_profiles[c,] == 0) for c in range(0, cells.d.shape[0])])


class NewtonRhapsonEstimator:
    MAX_Q: float = 0.999

    def __init__(self, cells: CellsData, tree: EventTreeWithCounts):
        self.cells = cells
        self.tree = tree
        self.CN = create_cn_matrix(cells, tree)
        self.S = self.__create_S_matrix()
        self.S_0 = self.__create_S_0_matrix()
        self.SCN = self.CN * self.S

    def solve(self, params: Parameters) -> Parameters:
        p_vec = np.array([params.e, params.m, params.q])

        for i in range(0, 30):
            target = self.calculate_target(params)
            jacobian = self.calculate_jacobian(params)
            if np.sum(self.S_0) == 0.0:
                inverse_sub_jacobian = np.linalg.inv(jacobian[1:3, 1:3])
                inverse_jacobian = np.zeros((3,3))
                inverse_jacobian[1:3, 1:3] = inverse_sub_jacobian
            else:
                inverse_jacobian = np.linalg.inv(self.calculate_jacobian(params))
            diff_vec = np.matmul(inverse_jacobian, target)

            for i in range(0, diff_vec.shape[0]):
                if p_vec[i] - diff_vec[i] < 0:
                    diff_vec[i] = 0.7 * p_vec[i]
            p_vec = p_vec - diff_vec
            params.e = p_vec[0]
            params.m = p_vec[1]
            params.q = p_vec[2]
            params.q = min(params.q, NewtonRhapsonEstimator.MAX_Q)
        return params

    def __create_S_matrix(self) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        for cell in range(0, self.cells.d.shape[0]):
            result[cell,] = self.cells.cell_cluster_sizes[cell]
        return result

    def __create_S_0_matrix(self) -> np.ndarray:
        S = self.__create_S_matrix()
        S[self.CN > 0] = 0
        return S

    def __create_alpha_matrix(self, params: Parameters) -> np.ndarray:
        b = self.CN.copy().astype(float)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell] * params.m
        b = b + self.__create_S_0_matrix() * params.e
        return b

    def __create_b_special_matrix(self, alpha: np.ndarray, params: Parameters,
                                  func: Callable[[float, float], float]) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        coef = math.exp(math.log(1 - params.q) - math.log(params.q))
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = func(coef * alpha[cell, bin], self.cells.d[cell, bin])
        return result

    def calculate_target(self, p: Parameters) -> np.ndarray:
        A = self.__create_alpha_matrix(p)
        B_inv = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum(
            [1 / (a + i) for i in range(0, int(d))]))
        S_0 = self.S_0
        SCN = self.SCN

        log_coef = math.log(1 - p.q) - math.log(p.q)
        f_1 = 0.0
        if np.sum(S_0) > 0.0:
            f_1 = math.log(1 - p.q) * math.exp(log_coef + math.log(np.sum(S_0))) + \
              (0.0 if np.sum(S_0 * B_inv) == 0.0 else math.exp(log_coef + math.log(np.sum(S_0 * B_inv))))
        f_2 = math.log(1 - p.q) * math.exp(log_coef + math.log(np.sum(SCN))) + \
              math.exp(log_coef + math.log(np.sum(SCN * B_inv)))
        f_3 = -math.log(1 - p.q) * math.exp(math.log(np.sum(A)) - 2.0 * math.log(p.q)) \
              - math.exp(math.log(np.sum(A)) - math.log(p.q)) + \
              + math.exp(math.log(np.sum(self.cells.d)) - math.log(p.q)) \
              - math.exp(math.log(np.sum(A * B_inv)) - 2.0 * math.log(p.q))

        return np.array([f_1, f_2, f_3])

    def calculate_jacobian(self, p: Parameters) -> np.ndarray:
        A = self.__create_alpha_matrix(p)
        S_0 = self.S_0
        SCN = self.SCN
        B_inv = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum(
            [1 / (a + i) for i in range(0, int(d))]))
        B_inv_2 = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum(
            [1 / (a + i) ** 2 for i in range(0, int(d))]))

        log_coef = math.log(1 - p.q) - math.log(p.q)

        f_1_e = 0.0 if np.sum(S_0 * S_0 * B_inv_2) == 0.0 else -math.exp(
            math.log(np.sum(S_0 * S_0 * B_inv_2)) + 2.0 * log_coef)
        f_2_m = -math.exp(math.log(np.sum(SCN * SCN * B_inv_2)) + 2.0 * log_coef)

        f_1_q = 0.0
        if np.sum(S_0) > 0:
            f_1_q = -math.log(1 - p.q) * math.exp(math.log(np.sum(S_0)) - 2.0 * math.log(p.q)) \
                - math.exp(math.log(np.sum(S_0)) - math.log(p.q)) + \
                + 0.0 if np.sum(S_0 * A * B_inv_2) == 0.0 else math.exp(
            math.log(np.sum(S_0 * A * B_inv_2)) - 2.0 * math.log(p.q) + log_coef) + \
                                                               - 0.0 if np.sum(S_0 * B_inv) == 0.0 else math.exp(
            math.log(np.sum(S_0 * B_inv)) - 2.0 * math.log(p.q))

        f_2_q = -math.log(1 - p.q) * math.exp(math.log(np.sum(SCN)) - 2.0 * math.log(p.q)) \
                - math.exp(math.log(np.sum(SCN)) - math.log(p.q)) + \
                + math.exp(math.log(np.sum(SCN * A * B_inv_2)) - 2.0 * math.log(p.q) + log_coef) + \
                - math.exp(math.log(np.sum(SCN * B_inv)) - 2.0 * math.log(p.q))

        f_3_q = 2.0 * math.log(1 - p.q) * math.exp(math.log(np.sum(A)) - 3.0 * math.log(p.q)) + \
                + math.exp(math.log(np.sum(A)) - 2.0 * math.log(p.q) - math.log(1 - p.q)) + \
                + math.exp(math.log(np.sum(A)) - 2.0 * math.log(p.q)) \
                - math.exp(math.log(np.sum(self.cells.d)) - 2.0 * math.log(p.q)) \
                - math.exp(math.log(np.sum(A * A * B_inv_2)) - 4.0 * math.log(p.q)) + \
                + 2.0 * math.exp(math.log(np.sum(A * B_inv)) - 3.0 * math.log(p.q))

        f_1_m = 0.0
        f_2_e = 0.0
        f_3_m = f_2_q
        f_3_e = f_1_q

        return np.array([[f_1_e, f_1_m, f_1_q], [f_2_e, f_2_m, f_2_q], [f_3_e, f_3_m, f_3_q]])
