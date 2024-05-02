#   mhodlr

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Run MATLAB Script on GitHub-Hosted Runner](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml/badge.svg)](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml)

## Matrix computation in HODLR representation

This repository contains the code for HODLR matrix as well as its basic matrix computation.

Differential equations often result in rank-structured matrices associated with off-diagonal blocks of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.

We consider computation on the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix. HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. In this repository, we implement HODLR computations in MATLAB, which is aimed for convenient API for HODLR operations. 


Setup
-------

One can fork this repository, and simply download this repository via
```bash
git clone https://github.com/<username>/mhodlr.git
```
and run the command below:
```bash
cd mhodlr
```


Then, one can try this simply example to verify its functionality:
```matlab
A = spdiags(ones(n, 1) * [2 8 -1],  -1:1, n, n); % generate test matrix
hA = hodlr(A); % Convert A to HODLR format
rA = recover(A); % Reconstruct hA into dense matrix
disp(norm(recover(full(rA - A)),2)); % Test error
```

Also, simulating adaptive precisions for each level of compression is allowed, please use the @mphodlr for low precision HODLR simulation. Users may compare the following code for usage guidance:
```matlab
rng(0);
A = rand(15,15); % Generate 15 by 15 random matrix

% Create precisions for each level; Level 1 use precision u1, level 2 use precision u2, ...
u1 = precision('q52');
u2 = precision('q52');
u3 = precision('h');
u_chain = prec_chain(u1, u2, u3);

% Usual call for full working precision 
hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
rA = recover(hA)
norm(full(rA - A),2)  % Compute the error

% Call mixed precision HODLR representation
mphA = mphodlr(u_chain, A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
mprA = recover(hA)
norm(full(mprA - A),2) % Compute the error

```

We refer to [document](https://github.com/chenxinye/mhodlr/blob/main/docs/source/start.rst) for usage in detail.


Support functions
---------------

|  Matrix computations | API|
|  ----  | ----  |
| Matrix transpose   | [``hodlr().transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) [``mphodlr().transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m)|
| Matrix inversion   | [``hodlr().inverse()``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) [``mphodlr().inverse()``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m)|
| Matrix (vector) multiplication | [``hdot``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hdot.m) |
| LU factorization   | [``hlu``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hlu.m)|
| Cholesky factorization  | [``hchol``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hchol.m)|

Acknowledgement
---------------
This repository is maintained by [InEXASCALE](https://www.karlin.mff.cuni.cz/~carson/inexascale) team at Charles University. 
