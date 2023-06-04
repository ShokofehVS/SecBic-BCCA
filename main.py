import numpy as np
from BCCA import BiCorrelationClusteringAlgorithm
from data_processing import load_yeast_tavazoie
from SecBic_BCCA import SecuredBiCorrelationClusteringAlgorithm

# load yeast data used in the original Cheng and Church's paper
data = load_yeast_tavazoie().values

# missing value imputation suggested by Cheng and Church
missing = np.where(data < 0.0)
data[missing] = np.random.randint(low=0, high=800, size=len(missing[0]))



print(data[0:2])
print(data.shape)


print("secure BCCA:")

sec_bcca = SecuredBiCorrelationClusteringAlgorithm(correlation_threshold=0.9, min_cols=3, data=data[0:2])
sec_bcca.run()
