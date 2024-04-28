#   mhodlr

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Matrix computation in HODLR representation

This repository contains the code for HODLR matrix as well as its basic matrix computation.

Differential equations often result in rank-structured matrices, most of whose off-diagonal blocks are often of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.

We consider computation on the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix. HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. In this repository, we implement HODLR computations in MATLAB, which is aimed for convenient API for HODLR operations. 

-------------

Setup
-------

Simply download this repo, and run the command below
```bash
cd mhodlr
```

Support functions
---------------


Acknowledgement
---------------
This repository is maintained by InEXASCALE team at Charles University. 
