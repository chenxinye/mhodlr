Get stated
======================================

Simply construct HODLR representation for array A via 

.. code:: matlab

   A = spdiags(ones(n, 1) * [1 5  -1],  -1:1, 100, 100);
   hA = hodlr(A, 'qr', 1.0e-8, 2, 3); % or simply use ``hA = hodlr(A)`` by leaving other parameters as default

.. admonition:: Note

    You can speficy four paramters for your initialization.  The first parameter is the matrix to be converted, the second is to specify the method low rank approximation, which default is 'svd', the third parameter is the tolerance for the low rank approximation, the fourth and the fifth refer to the minimum block size and maximum level for the cluster tree, respectively.

