import math
import numpy as np
from Pyfhel import Pyfhel

from _base import BaseBiclusteringAlgorithm
from model import Bicluster, Biclustering
from itertools import combinations
from sklearn.utils.validation import check_array


class SecuredBiCorrelationClusteringAlgorithm(BaseBiclusteringAlgorithm):
    """Secure Bi-Correlation Clustering Algorithm (SecBCCA) using Pyfhel and CKKS scheme"""
    def __init__(self, correlation_threshold=0.9, min_cols=3, data=np.array([])):
        self.correlation_threshold = correlation_threshold
        self.min_cols = min_cols
        self.data = data 

    def run(self):
        """Compute biclustering.

        Parameters
        ----------
        
        """

        # Creating empty Pyfhel object
        HE = Pyfhel()
        ckks_params = {
            'scheme': 'CKKS',
            'n': 2 ** 15,
            'scale': 2 ** 30,
            'qi_sizes': [60] + 15 * [30] + [60]
        }
        HE.contextGen(**ckks_params)  # Generate context for ckks scheme
        HE.keyGen()  # Key Generation: generates a pair of public/secret keys
        HE.rotateKeyGen()  # Rotation values in the vector
        HE.relinKeyGen()  # Relinearization key generation

        data = check_array(self.data, dtype=np.double, copy=True)
        self._validate_parameters()

        self.num_rows, self.num_cols = data.shape  # Moved inside the run method
        print(data.shape)

        # c_data = np.empty((self.num_rows, self.num_cols), dtype=object)
        # for i, row in enumerate(data):
        #     encrypted_row = [HE.encrypt(element) for element in row]
        #     c_data[i] = encrypted_row

        # # Encrypt the matrix
        # c_data = np.empty((self.num_rows, self.num_cols), dtype=object)
        # for i in range(data.shape[0]):
        #     for j in range(data.shape[1]):
        #         c_data[i, j] = HE.encrypt(data[i, j])
        # c_data = np.array(c_data, dtype=object)
        c_data = np.array([HE.encrypt(data[i]) for i in range(self.num_rows)])        
        biclusters = []
        print(c_data.shape)

        #########

        for i, j in combinations(range(self.num_rows), 2):
            cols, corr = self._encrypted_find_cols(HE, c_data[i], i, c_data[j], j, self.num_cols)
            print(corr)

            if len(cols) >= self.min_cols and corr >= self.correlation_threshold:
                rows = [i, j]

                for k, r in enumerate(c_data):
                    if k != i and k != j and self._accept(HE, rows, cols, r):
                        rows.append(k)

                b = Bicluster(rows, cols)

                if not self._exists(biclusters, b):
                    biclusters.append(b)

        #########
        
        for b in biclusters:
            print(b)

        return Biclustering(biclusters)
    
    def _encrypted_find_cols(self,HE, ri, row_i, rj, row_j, l_cols):
        """Finds the column subset for which the correlation between ri and rj
        stands above the correlation threshold.
        """
        cols = np.arange(l_cols, dtype=int)
        corr = self.pearson_corr_ckks(HE, ri, rj)

        while corr < self.correlation_threshold and len(cols) >= self.min_cols:
            imax = self._find_max_decrease(HE, ri, row_i, rj, row_j, cols)
            cols = np.delete(cols, imax)
            n_ri = HE.encrypt(self.data[row_i][cols])
            n_rj = HE.encrypt(self.data[row_j][cols])
            corr = self.pearson_corr_ckks(HE, n_ri, n_rj)
        return cols, corr
    
    def _find_max_decrease(self,HE, ri, row_i, rj, row_j, indices):
        """Finds the column which deletion causes the maximum increase in
        the correlation value between ri and rj
        """
        kmax, greater = -1, float('-inf')

        for k in range(len(indices)):
            ind = np.concatenate((indices[:k], indices[k+1:]))
            n_ri = HE.encrypt(self.data[row_i][ind])
            n_rj = HE.encrypt(self.data[row_j][ind])
            result = self.pearson_corr_ckks(HE, n_ri, n_rj)

            if result > greater:
                kmax, greater = k, result

        return kmax
    
    def _accept(self, HE, rows, cols, r):
        """Checks if row r satisfies the correlation threshold."""
        for i in rows:
            corr = self.pearson_corr_ckks(HE, r, HE.encrypt(self.data[i, cols]))

            if corr < self.correlation_threshold:
                return False

        return True
    
    def _exists(self, biclusters, bic):
        """Checks if a bicluster has already been found."""
        for b in biclusters:
            if len(b.rows) == len(bic.rows) and len(b.cols) == len(bic.cols) and \
                 np.all(b.rows == bic.rows) and np.all(b.cols == bic.cols):
                return True
        return False

    
    def _pearson_corr_ckks(self, HE, ctxt_x, ctxt_y):
        # Compute the number of values in the input ciphertext arrays
        num_values_x = len(ctxt_x)
        num_values_y = len(ctxt_y)
        print("len of x: ", num_values_x)
        print("len of y: ", num_values_y)

        # Check if input arrays are empty
        if num_values_x == 0 or num_values_y == 0:
            raise ValueError("Input arrays are empty.")

        # Compute the mean value of ctxt_x
        ctxt_mean_x = HE.cumul_add(ctxt_x, in_new_ctxt=True) / num_values_x
        print("mean encrypted: ", HE.decrypt(ctxt_mean_x)[0])

        # Compute the mean value of ctxt_y
        ctxt_mean_y = HE.cumul_add(ctxt_y, num_values_y)/ num_values_y

        # Compute the covariance of ctxt_x and ctxt_y
        ctxt_x_cov = ctxt_x - ctxt_mean_x
        ctxt_y_cov = ctxt_y - ctxt_mean_y

        # Compute the element-wise product of ctxt_x_cov and ctxt_y_cov
        ctxt_xy_cov = ctxt_x_cov * ctxt_y_cov

        # Compute the sum of ctxt_xy_cov
        ctxt_xy_cov_sum = HE.cumul_add(ctxt_xy_cov, in_new_ctxt=True)

        # Compute the squared differences of ctxt_x_cov
        ctxt_x_sdiv = ctxt_x_cov ** 2

        # Compute the squared differences of ctxt_y_cov
        ctxt_y_sdiv = ctxt_y_cov ** 2

        # Compute the product of the sum of ctxt_x_sdiv and ctxt_y_sdiv
        ctxt_xy_sdiv = HE.cumul_add(ctxt_x_sdiv, in_new_ctxt=True) * HE.cumul_add(ctxt_y_sdiv, in_new_ctxt=True)

        # Compute the square root of ctxt_xy_sdiv
        sq_root = np.sqrt(HE.decrypt(ctxt_xy_sdiv)[0])
        print(sq_root)
        # Compute the final result by dividing ctxt_xy_cov_sum by ctxt_sq_root
        c_result = ctxt_xy_cov_sum / sq_root

        return c_result

    def sum_with_error(slef, vector, n):
        """ This function sums up a pyfhel objects that contains error values.

        :param vector: vector with error values that should be summed up
        :param n: number for correct values in the vector
        :return: Pyfehl object which contains the sum at the first position
        """

        result_vector = vector.copy()
        for i in range(n-1):
            vector_shift = vector << (i + 1)
            result_vector = result_vector + vector_shift

        return result_vector
    
    def pearson_corr_ckks(self, HE, ctxt_x, ctxt_y):

        """Calculates the Pearson correlation"""
        
        """1. x_cov"""
        ctxt_mean_x = HE.cumul_add(ctxt_x, in_new_ctxt=True) / self.num_cols
        # ctxt_sum_x = HE.encrypt(0)  # Initialize a ciphertext for the sum

        # for ctxt_elem_x in ctxt_x:
        #     ctxt_sum_x += ctxt_elem_x  # Add each element to the sum ciphertext

        # ctxt_mean_x = ctxt_sum_x / self.num_cols
        ctxt_x_cov = ctxt_x - ctxt_mean_x

        """2. y_cov"""

        ctxt_mean_y = HE.cumul_add(ctxt_y, in_new_ctxt=True) / self.num_cols
        # ctxt_sum_y = HE.encrypt(0)  # Initialize a ciphertext for the sum

        # for ctxt_elem_y in ctxt_y:
        #     ctxt_sum_y += ctxt_elem_y  # Add each element to the sum ciphertext

        # ctxt_mean_y = ctxt_sum_y / self.num_cols
        ctxt_y_cov = ctxt_y - ctxt_mean_y

        """3. xy_cov"""

        ctxt_xy_cov = ctxt_x_cov * ctxt_y_cov

        # relinearizate and rescale after ciphertext-ciphertext multiplication to reduce size and scale
        HE.rescale_to_next(ctxt_xy_cov)
        HE.relinKeyGen()
        ctxt_xy_cov = ~ctxt_xy_cov


        ctxt_xy_cov_sum = self.sum_with_error(ctxt_xy_cov, self.num_cols)


        """4/5. x_sdiv and y_sdiv"""
        ctxt_x_sdiv = ctxt_x_cov ** 2
        ctxt_y_sdiv = ctxt_y_cov ** 2

        # rescaling of results to reduce scale
        HE.rescale_to_next(ctxt_x_sdiv)
        HE.rescale_to_next(ctxt_y_sdiv)


        # relinearization of results to reduce size
        HE.relinKeyGen()
        ctxt_x_sdiv = ~ctxt_x_sdiv

        HE.relinKeyGen()
        ctxt_y_sdiv = ~ctxt_y_sdiv

        """6. xy_sdiv"""

        ctxt_xy_sdiv = self.sum_with_error(ctxt_x_sdiv, self.num_cols) * self.sum_with_error(ctxt_y_sdiv, self.num_cols)

        # relinearizate and rescale after ciphertext-ciphertext multiplication to reduce size and scale
        HE.rescale_to_next(ctxt_xy_sdiv)
        HE.relinKeyGen()
        ctxt_xy_sdiv = ~ctxt_xy_sdiv

        """7. sq_root"""

        # decrypt value
        ptxt_xy_sdiv = HE.decrypt(ctxt_xy_sdiv)[0]

        # calculate square root on plaintext value
        ptxt_sq_root = np.sqrt(ptxt_xy_sdiv)

        """8. r_xy"""

        # calculate results as ciphertext-plaintext division
        result_ctxt = ctxt_xy_cov_sum / ptxt_sq_root

        # decrypt & decode results
        result = np.round(HE.decrypt(result_ctxt)[0], decimals=6)

        return np.abs(result)

    

    def _validate_parameters(self):
        if not (0.0 <= self.correlation_threshold <= 1.0):
            raise ValueError("correlation_threshold must be >= 0.0 and <= 1.0, got {}".format(self.correlation_threshold))

        if self.min_cols < 3:
            raise ValueError("min_cols must be >= 3, got {}".format(self.min_cols))