#   mhodlr

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Matrix computation in HODLR representation

This repository contains the code for HODLR matrix as well as its basic matrix computation.

Differential equations often result in rank-structured matrices associated with off-diagonal blocks of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.

We consider computation on the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix. HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. In this repository, we implement HODLR computations in MATLAB, which is aimed for convenient API for HODLR operations. 


Setup
-------

Simply download this repository via
```bash
git clone https://github.com/chenxinye/mhodlr.git
```
and run the command below:
```bash
cd mhodlr
```


Then, one can try this simply example to verify its functionality:
```matlab
A = spdiags(ones(n, 1) * [2 8 -1],  -1:1, n, n); % generate test matrix
hA = hodlr(A); % Convert A to HODLR format
RA = recover(A); % Reconstruct hA into dense matrix
disp(norm(recover(full(RA - A)),2)); % Test error
```

We refer to [document](https://github.com/chenxinye/mhodlr/blob/main/docs/source/start.rst) for usage in detail.


Support functions
---------------


Acknowledgement
---------------
This repository is maintained by InEXASCALE team at Charles University. 
