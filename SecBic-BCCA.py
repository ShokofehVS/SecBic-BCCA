import math

import numpy as np
from Pyfhel import Pyfhel


def ckks_corr(v, w):
    """init of the CKKs scheme"""
    HE = Pyfhel()  # Creating empty Pyfhel object
    ckks_params = {
        'scheme': 'CKKS',  # can also be 'ckks'
        'n': 2 ** 14,  # Polynomial modulus degree. For CKKS, n/2 values can be
        #  encoded in a single ciphertext.
        #  Typ. 2^D for D in [10, 15]
        'scale': 2 ** 30,  # All the encodings will use it for float->fixed point
        #  conversion: x_fix = round(x_float * scale)
        #  You can use this as default scale or use a different
        #  scale on each operation (set in HE.encryptFrac)
        'qi_sizes': [60, 30, 30, 30, 60]  # Number of bits of each prime in the chain.
        # Intermediate values should be  close to log2(scale)
        # for each operation, to have small rounding errors.
    }
    HE.contextGen(**ckks_params)  # Generate context for ckks scheme
    HE.keyGen()  # Key Generation: generates a pair of public/secret keys
    HE.rotateKeyGen()

    """encode and encrypt v and w"""
    arr_v = np.array(v, dtype=np.float64)
    arr_w = np.array(w, dtype=np.float64)

    num_values_v = len(arr_v)
    num_values_w = len(arr_w)

    ctxt_v = np.array([HE.encrypt(arr_v[i]) for i in range(num_values_v)])
    ctxt_w = np.array([HE.encrypt(arr_w[i]) for i in range(num_values_w)])

    # ctxt_v = [[ 0.1  0.1  0.1 ...  0.1  0.1  0.1]
    #  [ 0.2  0.2  0.2 ...  0.2  0.2  0.2]
    #  [-0.3 -0.3 -0.3 ... -0.3 -0.3 -0.3]]

    """Calculates the Pearson correlation and returns its absolute value."""
    # vc = v - np.mean(v)

    mean_ctxt_v = sum(ctxt_v) / num_values_v

    # mean_ctxt_v = [-0.  0.  0. ...  0. -0. -0.]

    ctext_vc = ctxt_v - mean_ctxt_v

    # ctect_vc = [[ 0.1  0.1  0.1 ...  0.1  0.1  0.1][ 0.2  0.2  0.2 ...  0.2  0.2  0.2]
    # [-0.3 -0.3 -0.3 ... -0.3 -0.3 -0.3]]

    # wc = w - np.mean(w)
    mean_ctxt_w = sum(ctxt_w) / num_values_w

    # mean_ctxt_w = [1.83 1.83 1.83 ... 1.83 1.83 1.83]

    ctext_wc = ctxt_w - mean_ctxt_w

    # ctext_wc = [[-3.33 -3.33 -3.33 ... -3.33 -3.33 -3.33][ 0.47  0.47  0.47 ...  0.47  0.47  0.47]
    # [ 2.87  2.87  2.87 ...  2.87  2.87  2.87]]

    # data = np.array(np.round([HE.decrypt(ctext_wc[i]) for i in range(num_values_v)], decimals=2))

    # x = np.sum(vc * wc)
    mul_ctxt_wv = ctext_vc * ctext_wc

    # mul_ctxt_wv = [[-0.33 -0.33 -0.33 ... -0.33 -0.33 -0.33][ 0.09  0.09  0.09 ...  0.09  0.09  0.09]
    # [-0.86 -0.86 -0.86 ... -0.86 -0.86 -0.86]]

    np.array([HE.rescale_to_next(mul_ctxt_wv[i]) for i in range(num_values_w)])

    ctxt_x = sum(mul_ctxt_wv)

    # ctxt_x = [-1.1 -1.1 -1.1 ... -1.1 -1.1 -1.1]



    # y = np.sum(vc * vc) * np.sum(wc * wc)
    sq_ctxt_vc = ctext_vc ** 2
    sq_ctxt_wc = ctext_wc ** 2

    # sq_ctxt_vc = [[0.01 0.01 0.01 ... 0.01 0.01 0.01][0.04 0.04 0.04 ... 0.04 0.04 0.04][0.09 0.09 0.09 ... 0.09 0.09 0.09]]
    # sq_ctxt_wc = [[11.11 11.11 11.11 ... 11.11 11.11 11.11][ 0.22  0.22  0.22 ...  0.22  0.22  0.22][ 8.22  8.22  8.22 ...  8.22  8.22  8.22]]


    np.array([HE.rescale_to_next(sq_ctxt_vc[i]) for i in range(num_values_w)])
    np.array([HE.rescale_to_next(sq_ctxt_wc[i]) for i in range(num_values_w)])

    ctxt_y = sum(sq_ctxt_vc) * sum(sq_ctxt_wc)
    # ctxt_y = [2.77 2.74 2.74 ... 2.73 2.73 2.74] or [2.74 2.74 2.74 ... 2.61 2.74 2.75] or ...

    print(np.round(HE.decrypt(ctxt_y), decimals=2))

    HE.rescale_to_next(ctxt_y)

    # [-9.03  2.72  9.19 ...  2.84 56.98  4.09] !!!!!!!!!!! rescaling is not working ? is the size to big ?

    print(np.round(HE.decrypt(ctxt_y), decimals=2))


    # return np.abs(x / np.sqrt(y))
    # root_txt_y = ctxt_y ** 0.5
    root_txt_y = ctxt_y

    result_ctxt = ctxt_x - root_txt_y

    print(np.round(HE.decrypt(ctxt_x), decimals=2))
    print(np.round(HE.decrypt(root_txt_y), decimals=2))


    # TODO abs

    """Decrypt & Decode results"""

    result = np.round(HE.decrypt(result_ctxt), decimals=2)
    print(result)

    return result


h = np.array([0.1, 0.2, -0.3])
print(np.round(np.sum(h), decimals=3))
z = arr_y = np.array([-1.5, 2.3, 4.7])

print(ckks_corr(h, z))