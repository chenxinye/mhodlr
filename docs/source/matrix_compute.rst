Matrix operations
======================================

This section is concerned with the matrix computations with HODLR matrix. 

Matrix product
------------------------------------------------

Matrix-vector product and matrix-matrix product share the same rountine, one simply use ``hdot`` to perform comptation, the code example is as below:

.. code:: matlab

    rng(0);
    A = rand(100);
    x = rand(100, 1); % define vector
    X = rand(100, 80); % define matrix

    % Usual call for full working precision 
    hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
    rA = recover(hA);
    norm(rA - A, 2)  % Compute the error

    b = hdot(hA, x); 
    err = norm(recover(b) - A * x, 'fro');
    disp(err);

    b = hdot(hA, x, 'dense');
    err = norm(b - A * x, 'fro');
    disp(err);

    B = hdot(hA, X);
    err = norm(recover(B) - A * X, 'fro');
    disp(err);

    B = hdot(hA, X, 'dense');
    err = norm(B - A * X, 'fro');
    disp(err);

The third parameter is optional, which indicates whether or not the output is of hodlr format, one can also specify the parameter to `dense`. The holdr format sometimes requires to be receovered for further operation. 

The output:

.. code::

   1.4379e-13

   1.4379e-13

   2.6714e-12

   1.3396e-12





LU factorization
------------------------------------------------

Working precision
^^^^^^^^^^^^^^^^^^

The LU factorization is performed by the rountine ``routine``


.. code:: matlab

    % Output the factors L and U are hodlr format as default
    [L, U] = hlu(hA); 
    err = norm(hdot(L, U, 'dense') - A, 'fro');
    disp(err);

    % Set the factors L and U to the dense matrix format. 
    [L, U] = hlu(hA, 'dense');
    err = norm(L * U - A, 'fro');
    disp(err);


Same to the ``hdot``, the last parameter are used to specify whether or not the output are of hodlr format.

The output:

.. code:: 

   1.4124e-12

   1.3376e-12



Multiple precision
^^^^^^^^^^^^^^^^^^^^

The working preicion for LU factorization can be specified by the method ``mhlu``:

.. code:: matlab

    u = precision('h');
    [L, U] = mhlu(hA, u, 'hodlr');
    err = norm(hdot(L, U, 'dense') - A, 'fro');
    disp(err);


One can also load the mixed precision ``mhodlr`` objects via, for example:

.. code:: matlab

    u1 = precision('d');
    u2 = precision('s');
    u3 = precision('h');
    u4 = precision('b');

    u_chain = prec_chain(u1, u2, u3, u4);
    depth=5;
    eps=1e-5;
    aphA = amphodlr(u_chain, A, depth, 10, 'svd', eps); 
    mphA = mphodlr(u_chain, A, depth, 10, 'svd', eps); 

    u = precision('h'); % set the working precision to half
    [L, U] = mhlu(mphA, u, 'hodlr');
    err = norm(hdot(L, U, 'dense') - A, 'fro');
    disp(err);

    u = precision('s'); % set the working precision to single
    [L, U] = mhlu(aphA, u, 'hodlr');
    err = norm(hdot(L, U, 'dense') - A, 'fro');
    disp(err);