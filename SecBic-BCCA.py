import math
import numpy as np
from Pyfhel import Pyfhel


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

    # each entry of x and y is encrypted separately into a Pyfhel object
    ctxt_x = np.array([HE.encrypt(arr_x[i]) for i in range(num_values_x)])
    ctxt_y = np.array([HE.encrypt(arr_y[i]) for i in range(num_values_y)])

    """Calculates the Pearson correlation"""

    """1. x_cov"""

    ctxt_mean_x = sum(ctxt_x) / num_values_x
    ctxt_x_cov = ctxt_x - ctxt_mean_x

    """2. y_cov"""

    ctxt_mean_y = sum(ctxt_y) / num_values_y
    ctxt_y_cov = ctxt_y - ctxt_mean_y

    """3. xy_cov"""

    ctxt_xy_cov = ctxt_x_cov * ctxt_y_cov

    # relinearizate and rescale after ciphertext-ciphertext multiplication to reduce size and scale
    HE.relinKeyGen()
    ctxt_xy_cov = ~ctxt_xy_cov
    np.array([HE.rescale_to_next(ctxt_xy_cov[i]) for i in range(num_values_y)])

    ctxt_xy_cov_sum = sum(ctxt_xy_cov)

    """4/5. x_sdiv and y_sdiv"""
    ctxt_x_sdiv = ctxt_x_cov ** 2
    ctxt_y_sdiv = ctxt_y_cov ** 2

    # relinearization of results to reduce size
    HE.relinKeyGen()
    ctxt_x_sdiv = ~ctxt_x_sdiv

    HE.relinKeyGen()
    ctxt_y_sdiv = ~ctxt_y_sdiv

    # rescaling of results to reduce scale
    np.array([HE.rescale_to_next(ctxt_x_sdiv[i]) for i in range(num_values_y)])
    np.array([HE.rescale_to_next(ctxt_y_sdiv[i]) for i in range(num_values_y)])

    """6. xy_sdiv"""
    ctxt_xy_sdiv = sum(ctxt_x_sdiv) * sum(ctxt_y_sdiv)

    # relinearizate and rescale after ciphertext-ciphertext multiplication to reduce size and scale
    HE.relinKeyGen()
    ctxt_xy_sdiv = ~ctxt_xy_sdiv
    HE.rescale_to_next(ctxt_xy_sdiv)

    """7. sq_root"""

    # decrypt value
    ptxt_xy_sdiv = np.mean(HE.decrypt(ctxt_xy_sdiv))

    # calculate square root on plaintext value
    ptxt_sq_root = math.sqrt(ptxt_xy_sdiv)

    """8. r_xy"""

    # calculate results as ciphertext-plaintext division
    result_ctxt = ctxt_xy_cov_sum / ptxt_sq_root

    # decrypt & decode results
    result = np.round(np.mean(HE.decrypt(result_ctxt)), decimals=4)

    return result
