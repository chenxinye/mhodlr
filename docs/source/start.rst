Get Stated
======================================


Setup 
-------------

The support Matlab releases for running ``mhodlr``: MATLAB2023a, MATLAB2023b, MATLAB2024a, MATLAB2024b.

Platform Compatibility: Windows, macOS, Linux.


One can easily install ``mhodlr`` by simply downloading the software from the `GitHub Link (https://github.com/chenxinye/mhodlr) <https://github.com/chenxinye/mhodlr>`_

Or use the bash command if git is installed:

.. code:: bash

   git clone https://github.com/chenxinye/mhodlr.git


After downloading the package, in Matlab terminal, add enviromental variables like

.. code

   addpath('mhodlr/mhodlr')

Or one can put the code in the corresponding location and perform the computing. 

The alternative way to download is from `FileExchange <https://www.mathworks.com/matlabcentral/fileexchange/170891-mhodlr>`_.








Quick start
-------------



Simply construct HODLR representation for array A via 

.. code:: matlab

   A = spdiags(ones(10, 1) * [1 5  -1],  -1:1, 10, 10);
   hA = hodlr(A, 3, 5, 'svd', 1.0e-8); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

By printing out the variable ``hA``, one can obtain:

.. code:: matlab

   hA = 

   hodlr with properties:

                  U1: [5×1 double]
                  V2: [-1 0 0 0 0]
                  U2: [5×1 double]
                  V1: [0 0 0 0 1]
                  D: []
                  A11: [1×1 hodlr]
                  A22: [1×1 hodlr]
               level: 1
               shape: [10 10]
            max_level: 3
         bottom_level: 1
               vareps: 1.0000e-08
      min_block_size: 5



Then you can recover your HODLR representation by built-in function ``recover``. 

.. code:: matlab

   RA = recover(hA);
   err = norm(full(RA - A), 2); % check how much error of the approxiamtion
   disp(err);

.. admonition:: Note

    You can speficy four paramters for your initialization\: The first parameter is the matrix to be converted; the second and the third refer to the maximum level and minimum block size for the cluster tree, respectively; the fourth is to specify the method low rank approximation, which default is 'svd', the other options: 'qr'; the fifth parameter is the tolerance for the low rank approximation; the sixth parameter is used to determine the norm used for the truncation, where only two norms, i.e., Frobenius and 2 norm, are supported, and the Frobenius norm only used with ``method``='svd'. 

One can perform hierarchical matrix-vector product via ``hdot``

.. code:: matlab

   x = rand(10);
   b = hdot(A, x, 'dense');
   
   err = norm(b - x, 2); # Compute forward
   disp(err);
   
.. admonition:: Note

   The third parameter is specified as 'dense' to indicate the output is dense format, the other option is 'hodlr'.
