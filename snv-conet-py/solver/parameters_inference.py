import math
import random
from dataclasses import dataclass
from typing import Callable

import numpy as np

from generator import SNVModel
from solver import CellsData


@dataclass
class Parameters:
    m: float  # per allele coverage
    e: float  # sequencing error
    q: float  # read success probablity


def create_cn_matrix(cells: CellsData, model: SNVModel) -> np.ndarray:
    matrix = np.zeros(cells.d.shape)
    for cell in range(0, cells.d.shape[0]):
        matrix[cell, ] = model.node_to_cn_profile[cells.attachment[cell]]
    return matrix.astype(int)


class MMEstimator:

    @staticmethod
    def estimate(cells: CellsData, model: SNVModel) -> Parameters:
        CN_profiles = create_cn_matrix(cells, model)
        zero_copy_bins = MMEstimator.__count_bin_cell_pairs_with_zero_copy(cells, CN_profiles)
        cn_sum = MMEstimator.__sum_copy_numbers(cells, CN_profiles)

        if zero_copy_bins == 0:
            raise RuntimeError(f"Can't estimate sequencing error - there are no full deletions in the dataset!")
        per_allele_coverage = np.sum(cells.d[CN_profiles != 0]) / cn_sum
        sequencing_error = np.sum(cells.d[CN_profiles == 0]) / zero_copy_bins

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
        return sum([cells.cell_cluster_sizes[c] * np.sum(CN_profiles[c, ] == 0) for c in range(0, cells.d.shape[0])])


class NewtonRhapsonEstimator:

    def __init__(self, cells: CellsData, model: SNVModel):
        self.cells = cells
        self.model = model
        self.CN = create_cn_matrix(cells, model)
        self.S = self.__create_S_matrix()
        self.S_0 = self.__create_S_0_matrix()
        self.SCN = self.CN * self.S

    def solve(self, params: Parameters) -> Parameters:
        p_vec = np.array([params.e, params.m, params.q])
        print(f"Starting parameters: {params.e}, {params.m}, {params.q}")

        for i in range(0, 30):
            target = self.calculate_target(params)
            jacobian = self.calculate_jacobian(params)
            inverse_jacobian = np.linalg.inv(jacobian)
            diff_vec = np.matmul(inverse_jacobian, target)
            print(f"Iter {i}, target {target},\n jacobian: {jacobian},\n diff {np.matmul(inverse_jacobian, target)},\n inverse {inverse_jacobian}")

            for i in range(0, diff_vec.shape[0]):
                if p_vec[i] - diff_vec[i] < 0:
                    diff_vec[i] = 0.7 * p_vec[i]
            p_vec = p_vec - diff_vec
            params.e = p_vec[0]
            params.m = p_vec[1]
            params.q = p_vec[2]
            params.q = min(params.q, 0.999)
            print(f"Parameters: {params.e}, {params.m}, {params.q}")
        return params

    def __create_S_matrix(self) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        for cell in range(0, self.cells.d.shape[0]):
            result[cell, ] = self.cells.cell_cluster_sizes[cell]
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

    def __create_b_special_matrix(self, alpha: np.ndarray, params: Parameters, func: Callable[[float, float], float]) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        coef = math.exp(math.log(1 - params.q) - math.log(params.q))
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = func(coef * alpha[cell, bin], self.cells.d[cell, bin])
        return result

    def calculate_target(self, p: Parameters) -> np.ndarray:
        A = self.__create_alpha_matrix(p)
        B_inv = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum([1 / (a + i) for i in range(0, int(d))]))
        S_0 = self.S_0
        SCN = self.SCN

        log_coef = math.log(1 - p.q) - math.log(p.q)

        f_1 = math.log(1 - p.q) * math.exp(log_coef + math.log(np.sum(S_0))) + \
              math.exp(log_coef + math.log(np.sum(S_0 * B_inv)))
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
        B_inv = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum([1 / (a + i) for i in range(0, int(d))]))
        B_inv_2 = self.__create_b_special_matrix(A, p, lambda a, d: 0.0 if d == 0.0 or a == 0.0 else sum([1 / (a + i) ** 2 for i in range(0, int(d))]))

        log_coef = math.log(1 - p.q) - math.log(p.q)

        f_1_e = -math.exp(math.log(np.sum(S_0 * S_0 * B_inv_2)) + 2.0 * log_coef)
        f_2_m = -math.exp(math.log(np.sum(SCN * SCN * B_inv_2)) + 2.0 * log_coef)

        f_1_q = -math.log(1 - p.q) * math.exp(math.log(np.sum(S_0)) - 2.0 * math.log(p.q)) \
                - math.exp(math.log(np.sum(S_0)) - math.log(p.q)) + \
                + math.exp(math.log(np.sum(S_0 * A * B_inv_2)) - 2.0 * math.log(p.q) + log_coef) + \
                - math.exp(math.log(np.sum(S_0 * B_inv)) - 2.0 * math.log(p.q))

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



class QEstimator:
    def __init__(self, cells: CellsData, model: SNVModel, params: Parameters):
        self.cells = cells
        self.model = model
        self.cn_profiles = create_cn_matrix(cells, model)
        self.p = params
        self.means = self.create_mean_matrix()

    def solve(self, iters: int) -> float:
        q1 = math.log(self.p.q)
        best_q = (q1, self.get_log_lik(q1))
        for i in range(0, iters):
            q2 = q1 + random.uniform(-1.0, 1.0)
            if q2 >= 0.0:
                continue
            acc = self.get_acceptance_ratio(q1, q2)
            if math.log(random.uniform(0.0, 1.0)) < acc:
                log_lik = self.get_log_lik(q2)
                print(f"Accepted q = {q2} with log-lik: {log_lik} and acc {acc}")
                if best_q[1] < log_lik:
                    best_q = (q2, log_lik)
                q1 = q2
            if i % 10 == 0:
                print(f"Best: {best_q[0], best_q[1]}")

        return q1

    def get_acceptance_ratio(self, q1: float, q2: float) -> float:
        return self.get_log_lik(q2) - self.get_log_lik(q1)

    def get_log_lik(self, q2: float) -> float:
        means = np.sum(self.means)
        D = np.sum(self.cells.d)
        return np.sum(means) * (math.log(1 - math.exp(q2)) * math.exp(math.log(1 - math.exp(q2)) - q2)) \
               + np.sum(D) * q2 \
               + np.sum(self.create_special_matrix(q2))

    def create_mean_matrix(self) -> np.ndarray:
        def __create_S_matrix() -> np.ndarray:
            result = np.zeros(self.cells.d.shape)
            for cell in range(0, self.cells.d.shape[0]):
                result[cell,] = self.cells.cell_cluster_sizes[cell]
            return result

        S = __create_S_matrix()
        CN = self.cn_profiles
        result = S * CN * self.p.m
        result[CN == 0] = S[CN == 0] * self.p.e
        return result

    def create_special_matrix(self, q: float) -> np.ndarray:
        result = np.zeros(self.means.shape)
        log_coef = math.log(1 - math.exp(q)) - q

        def __special(a: float, d: int):
            if d == 0 or a == 0.0:
                return 0
            result = 0.0
            for i in range(0, int(d)):
                result += math.log(math.exp(log_coef + math.log(a)) + i)
            return result

        for i in range(0, result.shape[0]):
            for j in range(0, result.shape[1]):
                result[i, j] = __special(self.means[i, j], self.cells.d[i, j])
        return result