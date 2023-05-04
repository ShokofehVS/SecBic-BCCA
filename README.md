# SecBic-BCCA

**SecBic-BCCA**: **Sec**ured **Bic**lusterings - **B**i-**C**orrelation **C**lustering **A**lgorithm: privacy-preserving gene expression data analysis by biclustering algorithm -- bi-correlation clustering algorithm -- over gene expression data with Homomorphic Encryption operations such as sum, or matrix multiplication in Python under the MIT license.

We apply [Pyfhel](https://pyfhel.readthedocs.io/en/latest/) as a python wrapper for the Microsoft SEAL library on top of the existing implementation of the algorithm in [biclustlib](https://github.com/padilha/biclustlib/). 

## Installation
First you need to ensure that all packages have been installed.
+ See `requirements.txt`
+ numpy>=1.23.1
+ setuptools>=65.5.0
+ pandas>=1.5.0
+ scikit-learn>=1.1.1
+ Pyfhel>=3.3.1
+ matplotlib>=3.5.2
+ scipy>=1.9.0
+ munkres>=1.1.4

You can clone this repository:

	   > git clone https://github.com/ShokofehVS/SecBic-BCCA.git

If you miss something you can simply type:

	   > pip install -r requirements.txt

If you have all dependencies installed:

	   > pip3 install .

To install Pyfhel, on Linux,`gcc6` for Python (`3.5+`) should be installed. (more information regarding [installation of Pyfhel ](https://github.com/ibarrond/Pyfhel))

	   > apt install gcc 

## Biclustering Algorithm
Biclustering or simultaneous clustering of both genes and conditions as a new paradigm was introduced by [Cheng and Church's Algorithm (CCA)](https://www.researchgate.net/profile/George_Church/publication/2329589_Biclustering_of_Expression_Data/links/550c04030cf2063799394f5e.pdf). The concept of bicluster refers to a subset of
genes and a subset of conditions with a high similarity score, which measures the coherence of the genes and conditions in the bicluster. It also returns the list of biclusters for the given data set. 

## Gene Expression Data Set
Our input data is *yeast Saccharomyces cerevisiae cell cycle* taken from [Tavazoie et al. (1999)](https://pubmed.ncbi.nlm.nih.gov/10391217/) which was used in the orginal study by [Cheng and Church](https://www.researchgate.net/profile/George_Church/publication/2329589_Biclustering_of_Expression_Data/links/550c04030cf2063799394f5e.pdf);

## External Evaluation Measure
To measure the similarity of encrypted biclusters with non-encrypted version, we use Clustering Error (CE) as an external evaluation measure that was proposed by [Patrikainen and Meila (2006)](http://ieeexplore.ieee.org/abstract/document/1637417/);
