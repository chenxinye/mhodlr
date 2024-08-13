



<p align="center">
 <img src="https://github.com/chenxinye/mhodlr/blob/main/data/lg.png?raw=true" alt="drawing" width="1180"/>
</p>

 <h1 align="center">
  mhodlr: Matrix computations in HODLR representation
 
[![License](https://img.shields.io/badge/License-BSD_3--Clause-lightblue.svg)](https://opensource.org/licenses/BSD-3-Clause)
[![MATLAB](https://github.com/chenxinye/mhodlr/actions/workflows/ci.yml/badge.svg)](https://github.com/chenxinye/mhodlr/actions/workflows/ci.yml)
[![View mhodlr on File Exchange](https://www.mathworks.com/matlabcentral/images/matlab-file-exchange.svg)](https://www.mathworks.com/matlabcentral/fileexchange/170891-mhodlr)
[![Documentation Status](https://readthedocs.org/projects/mhodlr/badge/?version=latest)](https://mhodlr.readthedocs.io/en/latest/?badge=stable)
[![Open in MATLAB Online](https://www.mathworks.com/images/responsive/global/open-in-matlab-online.svg)](https://matlab.mathworks.com/open/github/v1?repo=chenxinye/mhodlr&file=mhodlr)



</h1>




## Abstract


Differential equations often result in rank-structured matrices associated with low-rank off-diagonal blocks. These matrices are often represented in a hierarchical format, and their operation often results in fast arithmetic, e.g., matrix-vector product.  The hierarchical matrix [2] is a class of dense rank-structured matrices with a hierarchical low-rank off diagonal block structure, which frequently arises from finite element discretization of an elliptic PDE, radial basis function interpolation, and boundary integral equations. Hierarchical Off-Diagonal Low-Rank (HODLR) matrix, as a typical hierarchical matrix, is formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices.

It is widely known that low precision can reduce data communication and be more energy- and storage-efficient. Regarding the IEEE standard for floating point, single precision arithmetic can be twice as fast as double precision on specific hardware and the half-precision arithmetic achieves 4 times speedup over double precision. Using the software mhodlr, one can know what precisions are required for the HODLR matrix construction by evaluating their reconstruction error and computations error by simulating various precisions. This repository is concerned with HODLR matrix construction as well as basic matrix computations with HOLDR matrices, which aims to provide a convenient API for HODLR operations and efficient simulations for mixed-precision and adaptive precision HODLR matrix computing [1]. Our low precision arithmetic is simulated in terms of [4].  


Our software mainly contains three modules

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

#### Note these routines work for ``@hodlr``, ``@mphodlr``, and ``@amphodlr`` modules, 

|  Matrix computations | API|
|  ----  | ----  |
| Matrix transpose   | [``hodlr(A).transpose()``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hodlr.m) [``mphodlr(A).transpose()``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40mphodlr/mphodlr.m)|
| Matrix iverse | [``inverse(H)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/inverse.m)|
| Matrix (vector) multiplication | [``hdot(A/H, H/B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hdot.m), [``mphdot``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40mphodlr/mphdot.m) |
| LU factorization   | [``hlu(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hlu.m), [``mhlu(H/A, prec)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/mhlu.m)|
| Cholesky factorization  | [``hchol(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hchol.m), [``mhchol(H/A, prec)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/mhchol.m)|
| Triangular solver (Lower triangular solver LX=B, Upper triangular solver XU=B) |[``htrsl(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/htrsl.m), [``htrsu(B, H)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/htrsu.m)|
| Linear solver (Ax = b) |[``lu_solve(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/lu_solve.m)|

Contributions
---------------

Any forms of contributions are welcomed. Our documents are still in progress; feel free to pull request and submit issues for suggestions. Before contributing code, we suggest to contact the maintainers. The contact information of maintainers can be found in  [MaintainerList](https://mhodlr.readthedocs.io/en/latest/teams.html).


Acknowledgement
---------------

This project is supported by the European Union (ERC, [InEXASCALE](https://www.karlin.mff.cuni.cz/~carson/inexascale), 101075632). Views and opinions expressed are those of the authors only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting
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
