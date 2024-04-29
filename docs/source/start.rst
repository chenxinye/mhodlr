Get stated
======================================

Simply construct HODLR representation for array A via 

.. code:: matlab

   A = spdiags(ones(10, 1) * [1 5  -1],  -1:1, 10, 10);
   hA = hodlr(A, 'svd', 1.0e-8, 2, 3); % or simply use ``hA = hodlr(A)`` by leaving other parameters as default

Then you can recover your representation by 

.. code:: matlab

   RA = recover(hA);
   norm(full(RA - A),2) % check how much error it arises


.. admonition:: Note

    You can speficy four paramters for your initialization.  The first parameter is the matrix to be converted, the second is to specify the method low rank approximation, which default is 'svd', the third parameter is the tolerance for the low rank approximation, the fourth and the fifth refer to the minimum block size and maximum level for the cluster tree, respectively.

