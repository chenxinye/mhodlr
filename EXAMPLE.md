
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
