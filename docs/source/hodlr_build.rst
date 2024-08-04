Build HODLR matrix
======================================



.. code:: matlab

   A = spdiags(ones(10, 1) * [1 5  -1],  -1:1, 10, 10);
   hA = hodlr(A, 3, 5, 'svd', 1.0e-8); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

The generators U, V are represented in sparse format, the level indicates the current level of the HODLR matrix. 