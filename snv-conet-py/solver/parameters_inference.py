import math
import random
from dataclasses import dataclass

import numpy as np
import pandas.core.computation.ops

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
        matrix[cell,] = model.node_to_cn_profile[cells.attachment[cell]]
    return matrix.astype(int)


class MMEstimator:

    def estimate(self, cells: CellsData, model: SNVModel) -> Parameters:
        CN_profiles = create_cn_matrix(cells, model)
        zero_copy_bins = self.__count_bin_cell_pairs_with_zero_copy(cells, model)
        cn_sum = self.__sum_copy_numbers(cells, CN_profiles)
        cn_squares = self.__get_sum_of_cn_squares(cells, CN_profiles)

        per_allele_coverage = np.sum(cells.d[CN_profiles != 0]) / cn_sum
        sequencing_error = np.sum(cells.d[CN_profiles == 0]) / zero_copy_bins

        zero_copy_bins2 = 0
        for cell in range(0, cells.d.shape[0]):
            for bin in range(0, cells.d.shape[1]):
                if CN_profiles[cell, bin] == 0:
                    zero_copy_bins2 += cells.cell_cluster_sizes[cell] ** 2
        read_success_probability = 1 - np.sum(cells.d) / (
                np.sum(cells.d ** 2) - (per_allele_coverage ** 2) * cn_squares - (
                sequencing_error ** 2) * zero_copy_bins2)
        return Parameters(
            m=per_allele_coverage,
            e=sequencing_error,
            q=read_success_probability
        )

    def __get_sum_of_cn_squares(self, cells: CellsData, cluster_cn_profiles: np.ndarray) -> int:
        res = 0
        for cell in range(0, cells.d.shape[0]):
            res += np.sum(cluster_cn_profiles[cell,] ** 2) * cells.cell_cluster_sizes[cell] ** 2
        return res

    def __sum_copy_numbers(self, cells: CellsData, cluster_cn_profiles: np.ndarray) -> int:
        res = 0
        for cell in range(0, cells.d.shape[0]):
            res += np.sum(cluster_cn_profiles[cell,]) * cells.cell_cluster_sizes[cell]
        return res

    def __count_bin_cell_pairs_with_zero_copy(self, cells: CellsData, model: SNVModel) -> int:
        res = 0
        for cell in range(0, cells.d.shape[0]):
            cn_profile = model.node_to_cn_profile[cells.attachment[cell]]
            res += cells.cell_cluster_sizes[cell] * np.sum(cn_profile == 0)
        return res


class NewtonRhapsonEstimator:

    def __init__(self, cells: CellsData, model: SNVModel):
        self.cells = cells
        self.model = model
        self.cn_profiles = create_cn_matrix(cells, model)

    def solve(self, params: Parameters) -> Parameters:
        p_vec = np.array([params.e, params.m, params.q])

        for i in range(0, 100):
            print(f"Parans before: {params.e}, {params.m}, {params.q}")
            print(f"Iter {i}, diff {np.matmul(self.calculate_jacobian(params), self.calculate_target(params))}")
            p_vec = p_vec - np.matmul(self.calculate_jacobian(params), self.calculate_target(params))
            params.e = p_vec[0]
            params.m = p_vec[1]
            params.q = p_vec[2]
            print(f"Parans: {params.e}, {params.m}, {params.q}")
        return params

    def __create_S_matrix(self) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        for cell in range(0, self.cells.d.shape[0]):
            result[cell,] = self.cells.cell_cluster_sizes[cell]
        return result

    def __create_S_0_matrix(self) -> np.ndarray:
        S = self.__create_S_matrix()
        S[self.cn_profiles > 0] = 0
        return S

    def __create_alpha_matrix(self, params: Parameters) -> np.ndarray:
        b = self.cn_profiles.copy().astype(float)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell] * params.m
        b = b + self.__create_S_0_matrix() * params.e
        return b

    def __create_b_inv_matrix(self, alpha: np.ndarray, params: Parameters) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        coef = math.exp(math.log(1 - params.q) - math.log(params.q))
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = self.__special(coef * alpha[cell, bin], self.cells.d[cell, bin])
        return result

    def __create_b_inv_2_matrix(self, alpha: np.ndarray, params: Parameters) -> np.ndarray:
        result = np.zeros(self.cells.d.shape)
        coef = math.exp(math.log(1 - params.q) - math.log(params.q))
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = self.__special2(coef * alpha[cell, bin], self.cells.d[cell, bin])
        return result

    def calculate_target(self, p: Parameters) -> np.ndarray:
        A = self.__create_alpha_matrix(p)
        S = self.__create_S_matrix()
        S_0 = self.__create_S_0_matrix()
        SCN = self.cn_profiles * S
        B_inv = self.__create_b_inv_matrix(A, p)

        f_1 = math.log(1 - p.q) * np.sum(S_0) + np.sum(S_0 * B_inv)
        f_2 = math.log(1 - p.q) * np.sum(SCN) + np.sum(SCN * B_inv)
        f_3 = (math.log(1 - p.q) + p.q) * np.sum(A) - p.q * np.sum(self.cells.d) + np.sum(A * B_inv)

        return np.array([f_1, f_2, f_3])

    def calculate_jacobian(self, p: Parameters) -> np.ndarray:
        A = self.__create_alpha_matrix(p)
        S = self.__create_S_matrix()
        S_0 = self.__create_S_0_matrix()
        SCN = self.cn_profiles * S
        B_inv = self.__create_b_inv_matrix(A, p)
        B_inv_2 = self.__create_b_inv_2_matrix(A, p)

        coef = math.exp(math.log(1 - p.q) - math.log(p.q))
        coef_inv = 1 / coef

        f_1_e = -coef * np.sum(S_0 * S_0 * B_inv_2)
        f_2_m = -coef * np.sum(SCN * SCN * B_inv_2)
        f_1_q = -1 / (1 - p.q) * np.sum(S_0) + 1 / p.q ** 2 * np.sum(S_0 * A * B_inv_2)
        f_2_q = -1 / (1 - p.q) * np.sum(SCN) + 1 / p.q ** 2 * np.sum(SCN * A * B_inv_2)
        f_3_q = -coef_inv * np.sum(A) - np.sum(self.cells.d) + 1 / p.q ** 2 * np.sum(A * A * B_inv_2)
        f_1_m = 0.0
        f_2_e = 0.0
        f_3_m = (math.log(1 - p.q) + p.q) * np.sum(SCN) + np.sum(SCN * B_inv) - coef * np.sum(SCN * A * B_inv_2)
        f_3_e = (math.log(1 - p.q) + p.q) * np.sum(S_0) + np.sum(S_0 * B_inv) - coef * np.sum(S_0 * A * B_inv_2)

        jacobian = np.array([[f_1_e, f_1_m, f_1_q], [f_2_e, f_2_m, f_2_q], [f_3_e, f_3_m, f_3_q]])
        # jacobian = jacobian / self.cells.d.size
        return np.linalg.inv(jacobian)

    def __special2(self, a: float, d: int):
        if d == 0 or a == 0.0:
            return 0
        result = 0.0
        for i in range(0, int(d)):
            result += 1 / (a + i) ** 2
        return result

    def __special(self, a: float, d: int):
        if d == 0 or a == 0.0:
            return 0
        result = 0.0
        for i in range(0, int(d)):
            result += 1 / (a + i)
        return result


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
        result[result == 0] = S[CN == 0] * self.p.e
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
