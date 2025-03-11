



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
[![Hits](https://hits.seeyoufarm.com/api/count/incr/badge.svg?url=https%3A%2F%2Fgithub.com%2Fchenxinye%2Fmhodlr&count_bg=%23C550DA&title_bg=%23555555&icon=&icon_color=%23E7E7E7&title=hits&edge_flat=false)](https://github.com/chenxinye/mhodlr/)
</h1>




## Abstract

Differential equations often give rise to rank-structured matrices characterized by low-rank off-diagonal blocks. These matrices can be conveniently represented in a hierarchical format, enabling efficient arithmetic operations such as fast matrix-vector multiplication. Among these, hierarchical matrices [2] form a class of dense rank-structured matrices with a hierarchical low-rank off-diagonal block structure. Such matrices frequently emerge in applications like the finite element discretization of elliptic partial differential equations (PDEs), radial basis function interpolation, and boundary integral equations. A prominent example of hierarchical matrices is the Hierarchical Off-Diagonal Low-Rank (HODLR) matrix, which is constructed by hierarchically partitioning the matrix using a binary cluster tree. At each level of the tree, all off-diagonal blocks are represented as low-rank matrices, enabling efficient storage and computation.

<p align="left">
 <img src="https://raw.githubusercontent.com/chenxinye/mhodlr/refs/heads/main/docs/fancy_flowchart.png" alt="drawing" width="1380"/>
</p>


Leveraging low-precision arithmetic has become increasingly important in computational mathematics due to its potential to reduce data communication, energy consumption, and storage requirements. According to the IEEE floating-point standard, single-precision arithmetic can achieve up to twice the speed of double precision on certain hardware, while half-precision arithmetic can offer up to a fourfold speedup. The trade-offs between precision and computational accuracy are critical when designing algorithms for high-performance computing.

The mhodlr software repository addresses this by providing tools to evaluate and optimize precision requirements for HODLR matrix construction. By simulating various levels of precision, the software allows users to assess reconstruction errors and computational errors associated with HODLR matrix operations. The repository also offers a user-friendly API for constructing HODLR matrices and performing basic matrix computations, facilitating simulations for mixed-precision and adaptive-precision arithmetic in HODLR matrix computing [1]. The low-precision arithmetic simulations employed in this work are based on the methodologies described in [4], ensuring both accuracy and efficiency.


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

Then add the destination folder where ``mhodlr`` is installed to MATLAB’s search path:
```bash
addpath('mhodlr/mhodlr')
```


# Examples
---------------

### Build HODLR matrices
```matlab
rng(0); %fix randomness
A = rand(500);
depth = 5;
min_block_size = 2;
epsilon = 1e-8;
hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
```

```matlab
% define the precisions
u1 = precision('d');
u2 = precision('s');
u3 = precision('h');
u4 = precision('b');
u5 = precision('q52');

% build the collection of precisions
u_chain = prec_chain(u1, u2, u3, u4, u5); % the order matters!

% build the precisions according to u_chain (each layer uses the corresponding preicison)
mphA = mphodlr(u_chain, A, depth, min_block_size, 'svd', epsilon);
mprA = recover(mphA); % recover from the HODLR format

% build the precisions automatically
aphA = amphodlr(u_chain, A, depth, min_block_size, 'svd', epsilon);
aprA = recover(aphA); % recover from the HODLR format
```

Simple example on usage is referred to  [EXAMPLE](https://github.com/chenxinye/mhodlr/blob/main/EXAMPLE.md).

For detailed of matrix computations, please check [docs](https://mhodlr.readthedocs.io/en/latest/matrix_compute.html) for detail. 

Basic support routines
---------------

#### Note these routines work for ``@hodlr``, ``@mphodlr``, and ``@amphodlr`` modules, 

|  Matrix computations | API|
|  ----  | ----  |
| Matrix summation  | [``add(A, B, ...)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/add.m)|
| Matrix substraction | [``sub(A, B, ...)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/sub.m)|
| Matrix transpose   | [``H.transpose()``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hodlr.m)|
| Matrix inverse | [``inverse(H)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/inverse.m)|
| Matrix (vector) multiplication | [``dot(A/H, H/B, ...)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/hdot.m), [``mdot(A/H, H/B, ...)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/mhdot.m) |
| LU factorization   | [``hlu(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hlu.m), [``mhlu(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/mhlu.m)|
| Cholesky factorization  | [``hchol(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hchol.m), [``mhchol(H/A)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/mhchol.m)|
| QR factorization  | [``hqr(H, method)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/hqr.m), [``mhqr(H, method)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/mhqr.m)|
| Triangular solver (Lower triangular solver LX=B, Upper triangular solver XU=B) |[``htrsl(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/htrsl.m), [``htrsu(B, H)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/htrsu.m)|
| Linear solver (Ax = b) |[``lu_solve(H, B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/lu_solve.m) [``qr_solve(method, H, B)``](https://github.com/chenxinye/mhodlr/blob/main/mhodlr/%40hodlr/qr_solve.m)|

Multiple precision routine
---------------
``mhodlr`` enables multiple precision control for matrix computation based on HODLR representation. It allows users to control the precision in a global environment, without the need to specify the precision everytime when the function is called. 

The precision setting is performed by 
```matlab
u = precision('h'); % or ther precision customization
set_prec(u);
```

Then all the HODLR arithmetic function starting with ''m'' will performed in precision ``u``.



Contributions
---------------

Any forms of contributions are welcomed. Our documents are still in progress; feel free to pull request and submit issues for suggestions. Before contributing code, we suggest to contact the maintainers. The contact information of maintainers can be found in  [MaintainerList](https://mhodlr.readthedocs.io/en/latest/teams.html).


Acknowledgement
---------------


This project is supported by the European Union (ERC, [InEXASCALE](https://www.karlin.mff.cuni.cz/~carson/inexascale), 101075632). Views and opinions expressed are those of the authors only and do not necessarily reflect those of the European Union or the European Research Council. Neither the European Union nor the granting
authority can be held responsible for them.


[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.13335429.svg)](https://doi.org/10.5281/zenodo.13335429)


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
