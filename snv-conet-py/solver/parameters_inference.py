import math
from dataclasses import dataclass

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
            print(f"Iter {i}, diff { np.matmul(self.calculate_jacobian(params), self.calculate_target(params))}")
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

    def calculate_target(self, params: Parameters) -> np.ndarray:
        b = self.__create_beta_matrix(params)
        a = self.__create_alpha_matrix(params)
        S = self.__create_S_matrix()
        S_CN = self.cn_profiles * S
        special_b = self.__create_special(b, self.cells.d)
        special_a = self.__create_special(a, self.cells.d)

        f_1 = math.log(1 - params.q) * np.sum(S[self.cn_profiles == 0]) + np.sum(
            S[self.cn_profiles == 0] * special_b[self.cn_profiles == 0])
        f_2 = math.log(1 - params.q) * np.sum(S_CN) + np.sum(S_CN * special_a)
        f_3 = (math.log(1 - params.q) + params.q) * (params.m * np.sum(S_CN) + params.e * np.sum(self.cn_profiles == 0)) \
              - params.q * np.sum(self.cells.d) \
              + params.m * np.sum(S_CN ** special_a) \
              + params.e * np.sum(S[self.cn_profiles == 0] * special_b[self.cn_profiles == 0])
        return np.array([f_1, f_2, f_3]) / self.cells.d.size

    def calculate_jacobian(self, p: Parameters) -> np.ndarray:
        b = self.__create_beta_matrix(p)
        a = self.__create_alpha_matrix(p)
        S = self.__create_S_matrix()
        S_CN = self.cn_profiles * S
        special_b = self.__create_special(b, self.cells.d)
        special_a = self.__create_special(a, self.cells.d)
        special2_b = self.__create_special2(b, self.cells.d)
        special2_a = self.__create_special2(a, self.cells.d)

        coef = math.exp(math.log(1 - p.q) - math.log(p.q))
        coef_inv = 1 / coef

        f_1_e = -coef * np.sum(S[self.cn_profiles == 0] * S[self.cn_profiles == 0] * special2_b[self.cn_profiles == 0])
        f_2_m = -coef * np.sum(S_CN * S_CN * special2_a)
        f_1_q = -1 / (1 - p.q) * np.sum(S[self.cn_profiles == 0]) + 1 / p.q ** 2 * p.e * np.sum(
            S * S * special2_b)
        f_2_q = -1 / (1 - p.q) * np.sum(S_CN) + 1 / p.q ** 2 * p.m * np.sum(
            S_CN * S_CN * special2_a)

        f_3_q = -coef_inv * (np.sum(S_CN) * p.m + np.sum(S[self.cn_profiles == 0]) * p.e) \
                - np.sum(self.cells.d) \
                + 1 / p.q ** 2 * np.sum(S_CN * S_CN * special2_a) * p.m ** 2 \
                + 1 / p.q ** 2 * np.sum(S[self.cn_profiles == 0] * S[self.cn_profiles == 0] * special2_b[self.cn_profiles == 0]) * p.e ** 2
        f_1_m = 0.0
        f_2_e = 0.0
        f_3_m = (math.log(1 - p.q) + p.q) * np.sum(S_CN) + np.sum(S_CN * special_a) -coef * p.m * np.sum(S_CN * S_CN * special2_a)
        f_3_e = (math.log(1 - p.q) + p.q) * np.sum(S[self.cn_profiles == 0]) + np.sum(S[self.cn_profiles == 0] * special_b[self.cn_profiles == 0]) \
                    -coef * p.e * np.sum(S[self.cn_profiles == 0] * S[self.cn_profiles == 0] * special2_b[self.cn_profiles == 0])

        jacobian = np.array([[f_1_e, f_1_m, f_1_q], [f_2_e, f_2_m, f_2_q], [f_3_e, f_3_m, f_3_q]])
        jacobian = jacobian / self.cells.d.size
        return np.linalg.inv(jacobian)

    def __create_beta_matrix(self, params: Parameters) -> np.ndarray:
        c = math.log(params.e) + math.log(1 - params.q) - math.log(params.q)
        c = math.exp(c)
        b = np.full(self.cells.d.shape, c)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell]
        return b

    def __create_alpha_matrix(self, params: Parameters) -> np.ndarray:
        c = math.log(params.m) + math.log(1 - params.q) - math.log(params.q)
        c = math.exp(c)

        b = self.cn_profiles.copy().astype(float)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell]
        b *= c
        return b

    def __special(self, a: float, d: int):
        if d == 0 or a == 0.0:
            return 0
        result = 0.0
        for i in range(0, int(d)):
            result += 1 / (a + i) ** 2
        return result

    def __special2(self, a: float, d: int):
        if d == 0 or a == 0.0:
            return 0
        result = 0.0
        for i in range(0, int(d)):
            result += 1 / (a + i)
        return result

    def __create_special(self, b, d):
        result = np.zeros(d.shape)
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = self.__special(b[cell, bin], d[cell, bin])
        return result

    def __create_special2(self, b, d):
        result = np.zeros(d.shape)
        for cell in range(0, result.shape[0]):
            for bin in range(0, result.shape[1]):
                result[cell, bin] = self.__special2(b[cell, bin], d[cell, bin])
        return result