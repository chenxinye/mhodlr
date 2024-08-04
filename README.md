##   mhodlr: Matrix computations in HODLR representation

[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![Documentation Status](https://readthedocs.org/projects/mhodlr/badge/?version=latest)](https://mhodlr.readthedocs.io/en/latest/?badge=latest)
[![MATLAB](https://github.com/chenxinye/mhodlr/actions/workflows/ci.yml/badge.svg)](https://github.com/chenxinye/mhodlr/actions/workflows/ci.yml)

## Abstract

This repository contains simulation code for mixed-precision and adaptive precision HODLR matrix construction as well basic matrix computations with HOLDR matrices. 

Differential equations often result in rank-structured matrices associated with off-diagonal blocks of low rank. These matrices are often represented in hierarchical format, etc., whose operation often results in fast arithmetic, e.g., matrix-vector product.  A hierarchical matrix is a class of dense rank-structured matrices with a hierarchical low-rank structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, and boundary integral equations. Technically, a matrix with a low-rank off-diagonal structure can be represented in a hierarchical matrix format.

HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. This repository is concerned with Hierarchical Off-Diagonal Low-Rank (HODLR) matrices; we implement HODLR computations in MATLAB, and aim to provide a convenient API for HODLR operations. We also provide mixed precision simulation code for HODLR matrix computing.   

Our software mainly contains three components

|  Class | Description|
|  ----  | ----  |
|  [``@hodlr``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hodlr.m) | Compute HODLR matrix|
|  [``@mphodlr``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40mphodlr/mphodlr.m) | Compute HODLR matrix in mixed precision (precisions are defined by the users) |
|  [``@amphodlr``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40amphodlr/amphodlr.m) | Compute HODLR matrix in adaptive precision (precisions are provided by the users) |



Setup
-------

The environment for running ``mhodlr`` is MATLAB2023a, MATLAB2023b, MATLAB2024a, MATLAB2024b.

One can fork this repository, and simply download this repository via
```bash
git clone https://github.com/<username>/mhodlr.git
```
and run the command below:
```bash
cd mhodlr/mhodlr
```

Simple example on usage is referred to  [EXAMPLE](https://github.com/chenxinye/mhodlr/blob/main/EXAMPLE.md).

Support routines
---------------

#### Note these routines work for both hodlr, mphodlr, amphodlr class.

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
[1] C. Erin, X. Chen and X. Liu, Mixed precision HODLR matrices, arXiv:2407.21637, (2024), https://doi.org/10.48550/arXiv.2407.21637.

[2] S. B¨orm, L. Grasedyck, and W. Hackbusch, Introduction to hierarchical matrices with
 applications, Eng. Anal. Bound. Elem., 27 (2003), pp. 405–422, https://doi.org/10.1016/S0955-7997(02)00152-2.

[3] N. J. Higham and S. Pranesh, Simulating low precision floating-point arithmetic, SIAM J.
 Sci. Comput., 41 (2019), pp. C585–C602, https://doi.org/10.1137/19M1251308.


License
----------------

This project is licensed under the terms of the [![License](https://img.shields.io/badge/License-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause).
