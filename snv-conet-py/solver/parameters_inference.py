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
            print(f"Iter {i}")
            p_vec = p_vec - np.matmul(self.calculate_jacobian(params), self.calculate_target(params))
            params.e = p_vec[0]
            params.m = p_vec[1]
            params.q = p_vec[2]
        return params

    def calculate_target(self, params: Parameters) -> np.ndarray:
        b = self.__create_beta_matrix(params)
        a = self.__create_alpha_matrix(params)

        D_e = 0
        D_m = 0
        D_q = 0
        for cell in range(0, b.shape[0]):
            for bin in range(0, b.shape[1]):
                if self.cn_profiles[cell, bin] == 0:
                    D_e += self.cells.cell_cluster_sizes[cell] * (math.log(1-params.q) + self.__special2(b[cell, bin], self.cells.d[cell, bin]))
                    D_q +=   self.__special2(b[cell, bin], self.cells.d[cell, bin]) * self.cells.cell_cluster_sizes[cell] * params.e
                    D_q += self.cells.cell_cluster_sizes[cell] * params.e * (
                                math.log(1 - params.q) + params.q)

                D_m += self.cn_profiles[cell, bin] * self.cells.cell_cluster_sizes[cell] * (math.log(1-params.q) + self.__special2(a[cell, bin], self.cells.d[cell, bin]))
                D_q +=  self.__special2(a[cell, bin], self.cells.d[cell, bin]) * \
                       self.cells.cell_cluster_sizes[cell] * params.m * self.cn_profiles[cell, bin]
                D_q += -params.q * self.cells.d[cell, bin]
                D_q += self.cn_profiles[cell, bin] * self.cells.cell_cluster_sizes[cell] * params.m * (math.log(1-params.q) + params.q)
        return np.array([D_e, D_m, D_q])

    def calculate_jacobian(self, params: Parameters) -> np.ndarray:
        b = self.__create_beta_matrix(params)
        a = self.__create_alpha_matrix(params)

        coef = (1- params.q) / params.q
        coef_inv = 1 / coef

        D_e_e = 0
        D_m_m = 0
        D_q_q = 0
        D_q_m = 0
        D_q_e = 0

        for cell in range(0, b.shape[0]):
            for bin in range(0, b.shape[1]):
                if self.cn_profiles[cell, bin] == 0:
                    D_e_e += - coef * self.cells.cell_cluster_sizes[cell] ** 2 * self.__special(b[cell, bin], self.cells.d[cell, bin])
                    D_q_q += 1/ (params.q ** 2) * (params.e * self.cells.cell_cluster_sizes[cell])**2 * self.__special(b[cell, bin], self.cells.d[cell, bin])
                    D_q_q += coef_inv * (params.e * self.cells.cell_cluster_sizes[cell])

                D_m_m += - coef * (self.cn_profiles[cell, bin] * self.cells.cell_cluster_sizes[cell])** 2 * self.__special(a[cell, bin], self.cells.d[cell, bin])
                D_q_q += 1/ (params.q ** 2) * (self.cn_profiles[cell, bin] * params.m * self.cells.cell_cluster_sizes[cell])**2 * self.__special(a[cell, bin], self.cells.d[cell, bin])
                D_q_q += self.cells.d[cell, bin]
                D_q_q += coef_inv * (self.cn_profiles[cell, bin] * params.m * self.cells.cell_cluster_sizes[cell])

                D_q_m += (self.cn_profiles[cell, bin] * self.cells.cell_cluster_sizes[cell]) * (-1/ (1- params.q)) + 1 / params.q ** 2 *  (self.cn_profiles[cell, bin] * self.cells.cell_cluster_sizes[cell]) ** 2 * params.m * self.__special(a[cell, bin], self.cells.d[cell, bin])
                D_q_e += (self.cells.cell_cluster_sizes[cell]) * (-1/ (1- params.q)) + 1 / params.q ** 2 * (self.cells.cell_cluster_sizes[cell]) ** 2 * params.e * self.__special(b[cell, bin], self.cells.d[cell, bin])
        jacobian = np.array([[D_e_e, 0.0, D_q_e], [0.0, D_m_m, D_q_m], [D_q_e, D_q_m, D_q_q]])
        return np.linalg.inv(jacobian)

    def __create_beta_matrix(self, params: Parameters) -> np.ndarray:
        b = np.full(self.cells.d.shape, params.e * (1 - params.q) / params.q)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell]
        return b

    def __create_alpha_matrix(self, params: Parameters) -> np.ndarray:
        b = self.cn_profiles.copy().astype(float)
        for cell in range(0, self.cells.d.shape[0]):
            b[cell,] *= self.cells.cell_cluster_sizes[cell]
        b *= params.m * (1 - params.q) / params.q
        return b

    def __special(self, a: float, d: int):
        if d == 0:
            return 0
        result = 0.0
        for i in range(0, int(d)):
            result += 1/ (a + i) ** 2
        return result

    def __special2(self, a: float, d: int):
        if d == 0:
            return 0
        result = 0.0
        for i in range(0, d):
            result += 1/ (a + i)
        return result
