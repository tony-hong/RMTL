# RMTL
An R Library for Multi-task Learning

# Description 
This package provides an efficient implementation of regularized multi-task learning comprising 10 algorithms applicable for regression, classification, joint feature selection, task clustering, low-rank learning, sparse learning and network incorporation. All algorithms are implemented basd on the accelerated gradient descent method and feature a complexity of O(1/k^2). Sparse model structure is induced by the solving the proximal operator.

# Installation
```R
#it will take a while
R CMD check ./
R CMD build ./
R CMD INSTALL RMTL_1.0.tar.gz
```
# Details
Please check ["RMTL-manuel.pdf"](https://github.com/transbioZI/RMTL/blob/master/RMTL-manual.pdf) for more details.

# Contact
If you have any question, please contact: hank9cao@gmail.com
