#   mhodlr: Matrix computation in HODLR representation

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Run MATLAB Script on GitHub-Hosted Runner](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml/badge.svg)](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml)

## Abstract

This repository contains the code for HODLR matrix as well as its basic matrix computation. 

Differential equations often result in rank-structured matrices associated with off-diagonal blocks of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.

 HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. This repository is concerned with the computation on the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix; we implement HODLR computations in MATLAB, which is aimed for convenient API for HODLR operations. Besides, we provide mixed precision simulation code for HODLR matrix computing.   


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

Examples
-----------

After get software downloaded, one can try this simply example to verify its functionality:
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
norm(rA - A,2)  % Compute the error

% Call mixed precision HODLR representation
mphA = mphodlr(u_chain, A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
mprA = recover(hA)
norm(mprA - A,2) % Compute the error

```

Regarding matrix computation, currently ``mhodlr`` supports matrix inversion, matrix transpose, matrix dot product, LU factorization, and Cholesky factorization.  
The matrix dot product (``hdot`` and ``mphdot``) supports inputs of hodlr class, array, or their mixed type. We showcase the matrix dot product below:

```matlab
A = rand(15,15);

hA = hodlr(A, 3, 2, 'svd'); 

C_appr1 = hdot(hA, A, 'double'); % output array format
C_appr2 = hdot(hA, hA); % output HODLR format
C_true = A * A;

norm(C_appr1 - A,2)
norm(recover(C_appr2) - A,2)
```

Similarly, for mixed precision operation, each level of the computation will follow the array ``u_chain`` precision settings:

```matlab
mphA = mphodlr(u_chain, A, 3, 2, 'svd');

mp_C_appr1 = mphdot(mphA, A, 'double'); % output array format
mp_C_appr2 = mphdot(mphA, mphA); % output HODLR format

norm(mp_C_appr1 - C_true,2)
norm(recover(mp_C_appr2) - C_true,2) 
```



For customized precision, we refer to [precisions](https://github.com/chenxinye/mhodlr/blob/main/docs/source/precision.rst) for definition. 

Also, we refer to [document](https://github.com/chenxinye/mhodlr/blob/main/docs/source/start.rst) for usage in detail.


Support functions
---------------

|  Matrix computations | API|
|  ----  | ----  |
| Matrix transpose   | [``hodlr().transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) [``mphodlr().transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m)|
| Matrix inversion   | [``hodlr().inverse()``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) [``mphodlr().inverse()``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m)|
| Matrix (vector) multiplication | [``hdot``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hdot.m) [``mphdot``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphdot.m) |
| LU factorization   | [``hlu``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hlu.m)|
| Cholesky factorization  | [``hchol``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hchol.m)|


Contributions
---------------
Any forms of contributions are welcomed. Our documents is still under writing, feel free to pull request and submit issues for suggestions. Before contributing code, we suggest contact the maintainers. The contact information of maintainers can be found in  [MaintainerList](https://github.com/chenxinye/mhodlr/blob/main/maintainerList).


Acknowledgement
---------------
This repository is maintained by [InEXASCALE](https://www.karlin.mff.cuni.cz/~carson/inexascale) team at Charles University. 


License
----------------

This project is licensed under the terms of the [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause).
