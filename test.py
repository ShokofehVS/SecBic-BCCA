import numpy as np
from BCCA import BiCorrelationClusteringAlgorithm
from data_processing import load_yeast_tavazoie
from SecBic_BCCA import SecuredBiCorrelationClusteringAlgorithm
from Pyfhel import Pyfhel
import math


# # Creating empty Pyfhel object
# HE = Pyfhel()
# ckks_params = {
#     'scheme': 'CKKS',
#     'n': 2 ** 15,
#     'scale': 2 ** 30,
#     'qi_sizes': [60] + 15 * [30] + [60]
# }
# HE.contextGen(**ckks_params)  # Generate context for ckks scheme
# HE.keyGen()  # Key Generation: generates a pair of public/secret keys
# HE.rotateKeyGen()  # Rotation values in the vector
# HE.relinKeyGen()  # Relinearization key generation

# # data = [10,10,10,10,10,10,10]
# # n = len(data)
# x=9
# c_x = HE.encrypt(x) 
# print("cipher text:", c_x)



def sqrt_approximation(x, a, num_terms):
    sqrt_a = math.sqrt(a)
    approximation = sqrt_a
    

    for i in range(1, num_terms):
        term = (x - a)**i/a**(i-1)
        term /= (2*i)*sqrt_a
        approximation += term

    return approximation

print(sqrt_approximation(9, 5, 30))
# print("normal mean : ", np.mean(data))

# print("cipher data: ")
# c_sum=02*i)
# for x in c_data:
#     c_sum+=x
#     print(HE.decrypt(x)[0])
    


# c_mean = c_sum / n

# print("encrupted mean : ", HE.decrypt(c_mean)[0])
