##   mhodlr: Matrix computations in HODLR representation

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Run MATLAB Script on GitHub-Hosted Runner](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml/badge.svg)](https://github.com/chenxinye/mhodlr/actions/workflows/myscript.yml)

## Abstract

This repository contains code for HODLR matrix construction as well basic matrix computations with HOLDR matrices. 

Differential equations often result in rank-structured matrices associated with off-diagonal blocks of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, and boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.


<img src=docs/demo.png width=300 />

HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. This repository is concerned with Hierarchical Off-Diagonal Low-Rank (HODLR) matrices; we implement HODLR computations in MATLAB, and aim to provide a convenient API for HODLR operations. We also provide mixed precision simulation code for HODLR matrix computing.   

Our software mainly contains three components


|  Class | Description|
|  ----  | ----  |
|  [``@hodlr``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) | Compute HODLR matrix|
|  [``@mphodlr``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m) | Compute HODLR matrix in mixed precision (precisions are defined by the users) |
|  [``@amphodlr``](https://github.com/chenxinye/mhodlr/blob/main/%40amphodlr/amphodlr.m) | Compute HODLR matrix in adaptive precision (precisions are provided by the users) |



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

After the software has been downloaded, users may need to be familiar with the precision class if mixed precision routine is needed, we refer to [precisions](https://github.com/chenxinye/mhodlr/blob/main/docs/source/precision.rst) for definitions. 

One can try this simple example to verify its functionality:
```matlab
A = spdiags(ones(n, 1) * [2 8 -1],  -1:1, n, n); % generate test matrix
hA = hodlr(A); % Convert A to HODLR format
rA = recover(hA); % Reconstruct hA into dense matrix
disp(norm(full(rA - A), 2) / norm(full(A), 2)); % relative test error
```

Also, simulating adaptive precisions for each level of compression is allowed; please use the @mphodlr for low precision HODLR simulations. Users may consult the following code for usage guidance:
```matlab
rng(0);
A = rand(15,15); % Generate 15 by 15 random matrix

% Create precisions for each level; Level 1 use precision u1, level 2 use precision u2, ...
u1 = precision('q43');
u2 = precision('q52');
u3 = precision('h');
u_chain = prec_chain(u1, u2, u3);

% Usual call for full working precision 
hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
rA = recover(hA)
norm(rA - A,2)  % Compute the error

% Call mixed precision HODLR representation
mphA = mphodlr(u_chain, A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
mprA = recover(mphA)
norm(mprA - A,2) % Compute the error

```
It should noted that we use 3 precisions for the HODLR matrix computation since there are three levels of the cluster tree. If the number of precisions is less than the depth of the cluster tree, the rest of the levels will be performed in the working precision. One can compute in a uniform precision for each level via 

```matlab
u_chain = prec_chain(u1, u1, u1); % for a hierachical cluster tree with depth of 3
``` 

Regarding matrix computations, currently ``mhodlr`` supports matrix inversion, matrix transpose, matrix dot product, LU factorization, and Cholesky factorization.  
The matrix dot product (``hdot`` and ``mphdot``) supports inputs of hodlr class, array, or their mixed type. We showcase the matrix dot product below:

```matlab
A = rand(15,15);

hA = hodlr(A, 3, 2, 'svd'); 

C_appr1 = hdot(hA, A, 'dense'); % output array format
C_appr2 = hdot(hA, hA); % output HODLR format
C_true = A * A;

norm(C_appr1 - C_true,2)
norm(recover(C_appr2) - C_true,2)
```

Similarly, for mixed precision operations, each level of the computation will follow the array ``u_chain`` precision settings:

```matlab
mphA = mphodlr(u_chain, A, 3, 2, 'svd');

mp_C_appr1 = mphdot(mphA, A, 'dense'); % output array format
mp_C_appr2 = mphdot(mphA, mphA); % output HODLR format

norm(mp_C_appr1 - C_true,2)
norm(recover(mp_C_appr2) - C_true,2) 
```

For reduced memory used, we suggest use shallow cluster tree for transforming.


Also, we refer to [document](https://github.com/chenxinye/mhodlr/blob/main/docs/source/start.rst) for usage in detail.


To do
---------------
[1] Adaptive cross approximation for HODLR format

[2] QR computing with HODLR format

[3] Matrix function computing with HODLR format

[4] Singular value decomposition

Support routines
---------------

#### Note these routines work for both hodlr and mphodlr class.

|  Matrix computations | API|
|  ----  | ----  |
| Matrix transpose   | [``hodlr(A).transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hodlr.m) [``mphodlr(A).transpose()``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphodlr.m)|
| Matrix (vector) multiplication | [``hdot(A/H, H/B)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hdot.m) [``mphdot``](https://github.com/chenxinye/mhodlr/blob/main/%40mphodlr/mphdot.m) |
| LU factorization   | [``hlu(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hlu.m)|
| Cholesky factorization  | [``hchol(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/hchol.m)|
| Triangular solver (Lower triangular solver LX=B, Upper triangular solver XU=B) |[``htrsl(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/htrsl.m), [``htrsu(B, H)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/htrsu.m)|
| Linear solver (Ax = b) |[``lu_solve(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/%40hodlr/lu_solve.m)|

Contributions
---------------
Any forms of contributions are welcomed. Our documents are still in progress; feel free to pull request and submit issues for suggestions. Before contributing code, we suggest to contact the maintainers. The contact information of maintainers can be found in  [MaintainerList](https://github.com/chenxinye/mhodlr/blob/main/maintainerList).


Acknowledgement
---------------
This project is supported by the European Union (ERC, [InEXASCALE](https://www.karlin.mff.cuni.cz/~carson/inexascale), 101075632). Views and opinions
 expressed are those of the authors only and do not necessarily reflect those of the European
 Union or the European Research Council. Neither the European Union nor the granting
 authority can be held responsible for them.




References
---------------
[1] Carson E., Chen X., Liu X., 2024, arXiv, arXiv:2407.21637. doi:10.48550/arXiv.2407.21637

[2] S. B¨ orm, L. Grasedyck, and W. Hackbusch, Introduction to hierarchical matrices with
 applications, Eng. Anal. Bound. Elem., 27 (2003), pp. 405–422, https://doi.org/10.1016/S0955-7997(02)00152-2.

[3] N. J. Higham and S. Pranesh, Simulating low precision floating-point arithmetic, SIAM J.
 Sci. Comput., 41 (2019), pp. C585–C602, https://doi.org/10.1137/19M1251308.


License
----------------

This project is licensed under the terms of the [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause).
