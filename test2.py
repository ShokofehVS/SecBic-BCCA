import math
import numpy as np
from Pyfhel import Pyfhel



def sum_with_error(vector, n):
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


def pearson_corr_ckks(x, y):
    """init of the ckks scheme as described in the pyfhel documentation"""

    HE = Pyfhel()  # Creating empty Pyfhel object
    ckks_params = {
        'scheme': 'CKKS',  # setting scheme to CKKS
        'n': 2 ** 14,  # set Polynomial modulus degree
        'scale': 2 ** 30,  # use default scale
        'qi_sizes': [60, 30, 30, 30, 60]
    }

    HE.contextGen(**ckks_params)  # Generate context for ckks scheme
    HE.keyGen()  # Key Generation: generates a pair of public/secret keys
    HE.rotateKeyGen()

    """encoding and encryption of x and y"""

    # save x and y as np with float numbers
    arr_x = np.array(x, dtype=np.float64)
    arr_y = np.array(y, dtype=np.float64)

    # save length of arrays
    num_values_x = len(arr_x)
    num_values_y = len(arr_y)

    # encrypt array
    ctxt_x = HE.encrypt(arr_x)
    ctxt_y = HE.encrypt(arr_y)

    """Calculates the Pearson correlation"""

    """1. x_cov"""

    ctxt_mean_x = HE.cumul_add(ctxt_x, in_new_ctxt=True) / num_values_x
    ctxt_x_cov = ctxt_x - ctxt_mean_x

    """2. y_cov"""

    ctxt_mean_y = HE.cumul_add(ctxt_y, in_new_ctxt=True) / num_values_y
    ctxt_y_cov = ctxt_y - ctxt_mean_y

    """3. xy_cov"""

    ctxt_xy_cov = ctxt_x_cov * ctxt_y_cov

    # relinearizate and rescale after ciphertext-ciphertext multiplication to reduce size and scale
    HE.rescale_to_next(ctxt_xy_cov)
    HE.relinKeyGen()
    ctxt_xy_cov = ~ctxt_xy_cov


    ctxt_xy_cov_sum = sum_with_error(ctxt_xy_cov, num_values_x)


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

    ctxt_xy_sdiv = sum_with_error(ctxt_x_sdiv, num_values_x) * sum_with_error(ctxt_y_sdiv, num_values_y)

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
    result = np.round(HE.decrypt(result_ctxt)[0], decimals=4)

    return np.abs(result)


x=[30,40,6,4,3,20,100,1000]
y=[10,50,60,60,20,1,100,100]


print("pearson encrypted:",pearson_corr_ckks(x, y))

def _corr(v, w):
        """Calculates the Pearson correlation and returns its absolute value."""
        vc = v - np.mean(v)
        wc = w - np.mean(w)


        x = np.sum(vc * wc)
        y = np.sum(vc * vc) * np.sum(wc * wc)

        res = x/np.sqrt(y)

        return np.abs(res)

print("pearson normal: ",_corr(x, y))