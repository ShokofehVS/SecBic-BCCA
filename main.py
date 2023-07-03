import numpy as np
from BCCA import BiCorrelationClusteringAlgorithm
from data_processing import load_yeast_tavazoie
from SecBic_BCCA import SecuredBiCorrelationClusteringAlgorithm

# load yeast data used in the original Cheng and Church's paper
data = load_yeast_tavazoie().values

# missing value imputation suggested by Cheng and Church
missing = np.where(data < 0.0)
data[missing] = np.random.randint(low=0, high=800, size=len(missing[0]))


# Set the size of the dataset
num_genes = 100  # Number of genes
num_conditions = 17  # Number of conditions


# Generate the synthetic yeast-like dataset
dataset = np.random.randint(1, 800, size=(num_genes, num_conditions))

# Print the dataset shape
print("Dataset shape:", dataset.shape)

#print(dataset)

print("nomal BCCA:")

bcca = BiCorrelationClusteringAlgorithm(correlation_threshold=0.9, min_cols=3)
normal_biclusters = bcca.run(dataset)

# print("secure BCCA:")

# sec_bcca = SecuredBiCorrelationClusteringAlgorithm(correlation_threshold=0.9, min_cols=3, data=dataset)
# sec_biclusters = sec_bcca.run()


