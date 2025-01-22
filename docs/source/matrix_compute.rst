Matrix operations
======================================

This section is concerned with the matrix computations with HODLR matrix. 

.. admonition:: Note

    All the following, we recover the HODLR representation to dense matrix via ``dense(.)``, the alternative way is to use ``recover(.)`` or simply use the member function of the HODLR representation via ``H.dense`` or ``H.dense()``. 


Matrix transpose and inverse
------------------------------------------------

The transpose of the HODLR matrix can be simply performed by 

.. code:: matlab
    
    tA = hA.transpose();
    disp(norm(hadd(tA, A', '-', 'dense'), 'fro')); % print error



The inverse of the HODLR matrix is 

.. code:: matlab

    rng(0); %fix randomness
    A = rand(500);
    depth = 5;
    min_block_size = 2;
    epsilon = 1e-8;
    hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
    iA = inverse(hA);  % return an the inverse in hodlr format
    disp(norm(hdot(hA, iA, 'dense') - eye(500), 'fro'))  % print error


The second parameter will be defaulted as 'hodlr', which determine the output as hodlr format; One can set it to 'dense' to return a dense matrix, for example:

.. code:: matlab
    
    iA = inverse(hA, 'dense');  % return an the inverse in dense format
    disp(norm(iA * A - eye(500), 'fro')); % print error
    

We provide the optional parameters to specify the output format, namely, `hodlr` and `dense`, and the third parameter is the algorithm to perform.


Regarding the specified-precision inverse, one can use 

.. code:: matlab
    
    u = precision('s');
    iA = minverse(hA, u); 
    disp(norm(hdot(hA, iA, 'dense') - eye(50), 'fro'))

    iA = minverse(hA, u, 'dense');  % return an the inverse in dense format
    disp(norm(iA * A - eye(50), 'fro'))
    
which is similar to the native inverse, but the second parameter is specifed as the precision, while the other two remain the same. 

.. admonition:: Note
    
    Note that we provide explicit way to truncate HODLR matrix to rank-p HODLR matrix, to do that, one can simply perform 

    .. code:: matlab

        hA = htruncate(hA, 3) # return rank-3 HODLR matrix



.. admonition:: Note

    The HODLR matrix arithmetic for summation, subtraction, and multiplication can be completed via the shortcut ``add(A, B)``, ``sub(A, B)`` and ``dot(A, B)`` if one just want the HODLR output. If one want a dense output, please following the interfaces as below. 

Matrix summation and subtraction
------------------------------------------------

The summation and subtraction of HODLR matrices (A + B and A - B, A or B are not necessarily of hodlr format) are performed by the method ``hadd``. The are four dominant parameters, i.e., input matrix A, input matrix B, symbol ('+' for summation, and '-' for subtraction) denoting summation and subtraction, and output format ('dense' or 'hodlr', the 'hodlr' is default).
To perform the operation of summation (subtraction is similar), one can use

.. code:: 

    hadd(A, B, operation('+'/'-'), 'dense'/'hodlr')
    

For example: 

.. code:: matlab
    
    rng(0); %fix randomness
    A = rand(500);
    B = rand(500);

    depth = 5;
    min_block_size = 2;
    epsilon = 1e-8;
    hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default
    hB = hodlr(B, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

    disp(norm(hadd(hA, B, '-', 'dense') - (A-B), 'fro'));
    disp(norm(hadd(hA, hB, '-', 'dense') - (A-B), 'fro'));
    disp(norm(hadd(A, hB, '-', 'dense') - (A-B), 'fro'));
    disp(norm(dense(hadd(hA, B, '-', 'hodlr')) - (A-B), 'fro'));
    disp(norm(dense(hadd(hA, hB, '-', 'hodlr')) - (A-B), 'fro'));
    disp(norm(dense(hadd(A, hB, '-', 'hodlr')) - (A-B), 'fro'));

    disp(norm(hadd(hA, B, '+', 'dense') - (A+B), 'fro'));  % print error
    disp(norm(hadd(hA, hB, '+', 'dense') - (A+B), 'fro')); % print error
    disp(norm(hadd(A, hB, '+', 'dense') - (A+B), 'fro')); % print error
    disp(norm(dense(hadd(hA, B, '+', 'hodlr')) - (A+B), 'fro')); % print error
    disp(norm(dense(hadd(hA, hB, '+', 'hodlr')) - (A+B), 'fro')); % print error
    disp(norm(dense(hadd(A, hB, '+', 'hodlr')) - (A+B), 'fro')); % print error






Note that one should ensure the two inputs are same structure (e.g., same depth) if they both are HODLR format.

Matrix product
------------------------------------------------

Matrix-vector product and matrix-matrix product share the same rountine, one simply use ``hdot`` for working precision or ``mhdot`` for varying precision to perform comptation.

Working precision
^^^^^^^^^^^^^^^^^^

The code example for working precision is as below:

.. code:: matlab

    rng(0);
    A = rand(100);
    x = rand(100, 1); % define vector
    X = rand(100, 80); % define matrix

    % Usual call for full working precision 
    hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
    rA = dense(hA);
    disp(norm(rA - A, 2)); % print error

    b = hdot(hA, x); 
    err = norm(dense(b) - A * x, 'fro');
    disp(err); % print error
 
    b = hdot(hA, x, 'dense');
    err = norm(b - A * x, 'fro');
    disp(err); % print error

    B = hdot(hA, X);
    err = norm(dense(B) - A * X, 'fro');
    disp(err); % print error

    B = hdot(hA, X, 'dense');
    err = norm(B - A * X, 'fro');
    disp(err); % print error

The third parameter is optional, which indicates whether or not the output is of hodlr format, one can also specify the parameter to `dense`. The holdr format sometimes requires to be receovered for further operation. 



Multiple precision
^^^^^^^^^^^^^^^^^^^^

To simulate specific precision for matrix-matrix product or matrix-vector product, the above example can be simply modifed to: 


.. code:: matlab

    u = precision('h');

    rng(0);
    A = rand(100);
    x = rand(100, 1); % define vector
    X = rand(100, 80); % define matrix

    % Usual call for full working precision 
    hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
 
    b = mhdot(hA, x, u, 'dense');
    err = norm(b - A * x, 'fro');
    disp(err); % print error

    B = mhdot(hA, X, u);
    err = norm(dense(B) - A * X, 'fro');
    disp(err); % print error

    B = mhdot(hA, X, u, 'dense');
    err = norm(B - A * X, 'fro');
    disp(err); % print error




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
    disp(err); % print error


Same to the ``hdot``, the last parameter are used to specify whether or not the output are of hodlr format.

.. admonition:: Note

   Note that the factors L and U are block lower and upper triangular matrix. 


Multiple precision
^^^^^^^^^^^^^^^^^^^^

The working preicion for LU factorization can be specified by the method ``mhlu``:

.. code:: matlab

    u = precision('h');
    [L, U] = mhlu(hA, u, 'hodlr');
    err = norm(hdot(L, U, 'dense') - A, 'fro');
    disp(err); % print error


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
    disp(err); % print error




Cholesky factorization
------------------------------------------------

The Cholesky factorization can be used similarly to LU factorization, which is implemented by the method ``hchol``. The following example briefly illustrate the usage of ``hchol``.


Working precision
^^^^^^^^^^^^^^^^^^

.. code:: matlab

    rng(0);
    R = rand(100);
    A = R'*R; % Generate symmetric positive definite matrix

    % Usual call for full working precision 
    hA = hodlr(A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.

    R = hchol(hA); % return a block upper triangular HODLR matrix

    disp(norm(hdot(R.transpose(), R, 'dense') - A, 'fro')) % print error


The first and second input of ``hchol`` is the input HODLR matrix and format of output, respectively; The second input is optional, which is defaulted as ``hodlr`` if it is missing. 

To generate the dense output, simply use:

.. code:: matlab

    R = hchol(hA, 'dense'); % return a 
    dusp(norm(R'*R - A, 'fro')); % print error


Multiple precision
^^^^^^^^^^^^^^^^^^^^

The usage of ``mhchol`` is similar, it proceeds by simply adding one additional parameter to determine the user-specific working precision.

.. code:: matlab

    % Create precisions for each level; Level 1 use precision u1, level 2 use precision u2, ...
    u1 = precision('q43');
    u2 = precision('q52');
    u3 = precision('b');
    u4 = precision('s');
    u_chain = prec_chain(u1, u2, u3, u4);


    % Call mixed precision HODLR representation
    amphA = amphodlr(u_chain, A, 3, 2, 'svd'); % Use maxmium level of 3 and minimum block size of 2, and perform SVD (default) low rank approximation.
    amprA = dense(amphA);
    norm(amprA - A,2) % Compute the error

    R = mhchol(amphA, u4); % or R = mhchol(hA, u4);
    disp(norm(hdot(R.transpose(), R, 'dense') - A, 'fro')); % print error




Matrix QR factorization
------------------------------------------------

We provide three implementations for QR factorizations. 

.. code:: matlab

    [Q, R] = hqr(hA, 'lintner');
    [Q, R] = hqr(hA, 'bebendorf');
    [Q, R] = hqr(hA, 'kressner');
    


