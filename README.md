#   mhodlr

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)

## Matrix computation in HODLR representation

This repository contains the code for HODLR matrix as well as its basic matrix computation.

Differential equations often result in rank-structured matrices, most of whose off-diagonal blocks are often of low rank. These matrices are often represented in Block Low Rank (BLR) format, Hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation \cite{amda13}, boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format. In this paper, we consider computation on the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix in a mixed precision manner. HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. 

-------------

