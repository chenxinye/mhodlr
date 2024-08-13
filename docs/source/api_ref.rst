API Reference
======================================

@precision 
-----------

precision.m
^^^^^^^^^^^^^

.. code:: matlab

    Parameters
    --------------------
    base - array | str, default='h'
       For the string type, the arithmetic format supports:
       'q43', 'fp8-e4m3'       - NVIDIA quarter precision (4 exponent bits,
                                 3 significand (mantissa) bits),
       'q52', 'fp8-e5m2'       - NVIDIA quarter precision (5 exponent bits,
                                 2 significand bits),
       'b', 'bfloat16'         - bfloat16,
       'h', 'half', 'fp16'     - IEEE half precision (the default),
       's', 'single', 'fp32'   - IEEE single precision,
       'd', 'double', 'fp64'   - IEEE double precision,
       'c', 'custom'           - custom format.
      The custom (base 2) format is defined by array (2-element, [t, emax]), 
      where t is the number of bits in the significand
      (including the hidden bit) and emax is the maximum value of the
      exponent.  The minimu exponent is taken to be emin = 1 - emax and
      the IEEE floating-point number representation is assumed, so that
      emax and the number of  bits e in the exponent are related by
      emax = 2^(e-1) - 1. 

    round - int, default=1
       1: round to nearest using round to even last bit to break ties
          (the default);
       2: round towards plus infinity (round up);
       3: round towards minus infinity (round down);
       4: round towards zero;
       5: stochastic rounding - round to the next larger or next smaller
          f.p. (floating-point) number with probability proportional to
          1 minus the distance to those f.p. numbers;
       6: stochastic rounding - round to the next larger or next smaller 
          f.p. number with equal probability.

    subnormal - int
        specifies whether subnormal numbers are supported, if `subnormal=0`, 
        subnormals are flushed to zero:
            0 = do not support subnormals (the default for base='b', i.e., bfloat16 format),
            1 = support subnormals (the default if base is set to others).

    explim - int, default=1
        ``explim = 0`` make emax (the maximal exponent) for the specified arithmetic disabled, 
        so overflow, underflow, or subnormal numbers will be produced only if necessary 
        for the data type.  This parameter is for exploring
        low precisions independent of range limitations.

    flip - int, default=0
        Determine whether each element of the rounded value has a randomly chosen bit in 
        its significand flipped with a certain probability.

    prob - double, default=0.5
         if flip = 1 then each element of the rounded
        value has a randomly chosen bit in its significand flipped with probability ``prob``.

    randfunc - func, default=@(n) rand(n, 1)
        The random function for stochastic rounding. 
        If options.randfunc is supplied, then in stochastic rounding (modes
        5 and 6) the random numbers used for rounding will be generated
        using that function. It should be a function that has a single argument
        for the number of random numbers to generate and returns a vector of
        the random numbers. 
        

    Properties
    --------------------
    Same as parameters, except 

    u - double
        Unit roundoff computed for current floating point format.

@hodlr
-----------

hodlr.m
^^^^^^^^^^^^^

.. code:: matlab

    Parameters
    --------------------
    A - double | single
        Matrix to be converted.
        
    max_level - int, default=9999
        Maximum level for cluster tree.

    min_block_size - int, default=2
        The minimum size for HODLR blocks.

    method - str, default='svd'
        The method to perform compression for off-diagonal blocks.

    vareps - double, default=1.0e-12
        The vareps value used for truncation of low rank approximation.

    trun_norm_tp - str, default='2'
        Norm type for the the off-diagonal block truncation ``||A - B||_trun_norm_tp <= vareps * ||B||``.
        
        
    Properties
    --------------------
    U1, V2, U2, V1 - double 
        The right upper block matrix of each level, we have A12 = U1 * V2 and A21 = U2 * V1.

    A11, A22 - hodlr 
        The diagonal block matrix in HODLR format (access in the next level). 

    shape - array
        The shape of object in the current level.  

    level - int
        The level for cluster tree.
    
    max_level - int
        The maximum level of cluster tree after transformation.

@mphodlr
-----------


mphodlr.m
^^^^^^^^^^^

.. code:: matlab
    
    Parameters
    --------------------
    precs - cell
        The cell array that contains the precision used for compression of each level. 
        Each element is a precision class.

    A - double | single
        Matrix to be converted.
        
    max_level - int, default=9999
        Maximum level for cluster tree.

    min_block_size - int, default=2
        The minimum size for HODLR blocks.

    method - str, default='svd'
        The method to perform compression for off-diagonal blocks.

    vareps - double, default=1.0e-12
        The vareps value used for truncation of low rank approximation.

    trun_norm_tp - str, default='2'
        Norm type for the the off-diagonal block truncation ``||A - B||_trun_norm_tp <= vareps * ||B||``.
    
        
    Properties
    --------------------
    U1, V2, U2, V1 - double 
        The right upper block matrix of each level, we have A12 = U1 * V2 and A21 = U2 * V1.

    A11, A22 - hodlr 
        The diagonal block matrix in HODLR format (access in the next level). 

    shape - array
        The shape of object in the current level. 

    level - int
        The level for cluster tree.
    
    max_level - int
        The maximum level of cluster tree after transformation.




@amphodlr 
-----------

amphodlr.m
^^^^^^^^^^^

.. code:: matlab

    Parameters
    --------------------
    precs - cell
        The cell array that contains the precision used for compression of each level. 
        Each element is a precision class.

    A - double | single
        Matrix to be converted.
        
    max_level - int, default=9999
        Maximum level for cluster tree.

    min_block_size - int, default=2
        The minimum size for HODLR blocks.

    method - str, default='svd'
        The method to perform compression for off-diagonal blocks.

    vareps - double, default=1.0e-12
        The vareps value used for truncation of low rank approximation.

    trun_norm_tp - str, default='2'
        Norm type for the the off-diagonal block truncation ``||A - B||_trun_norm_tp <= vareps * ||B||``.
    
    Properties
    --------------------
    U1, V2, U2, V1 - double 
        The right upper block matrix of each level, we have A12 = U1 * V2 and A21 = U2 * V1.

    A11, A22 - hodlr 
        The diagonal block matrix in HODLR format (access in the next level). 

    shape - array
        The shape of object in the current level. 

    level - int
        The level for cluster tree.
    
    max_level - int
        The maximum level of cluster tree after transformation.



Compute rountines
-----------------------------

hadd.m
^^^^^^^^^^^

.. code:: matlab

    The function is used for the operation of summation or subtraction for HODLR matrix A and B.

    Parameters
    --------------------
    A - hodlr | double
        Input matrix - hodlr class / dense tyle.
  
    B - hodlr | double
        Input matrix - hodlr class / dense tyle.
    
    operator - str, default = '+'
        The operator of add ('+') or minus ('-'), string type.

    oformat - str, default = 'hodlr'
        The format of returns.
    
    Returns
    --------------------
    C - hodlr | double
        Return matrix in hodlr class or dense array.
  




hdot.m
^^^^^^^^^^^

.. code:: matlab

    Compute dot product of A and B.

    Parameters
    --------------------
    A - hodlr | double
        Matrix in HODLR format or double array.
    
    B - hodlr | double
        Matrix in HODLR format or double array.

    oformat - str, default='hodlr'
        Output format: 'hodlr' or 'dense'.
    


    Returns
    --------------------
    C - hodlr | double
        The matrix of product. 

mhdot.m
^^^^^^^^^^^^^

.. code:: matlab

    Compute dot product of A and B.

    Parameters
    --------------------
    A - hodlr | double
        Matrix in HODLR format or double array.
    
    B - hodlr | double
        Matrix in HODLR format or double array.

    prec - precision
        Precision for the matrix-vector product.

    oformat - str, default='hodlr'
        Output format: 'hodlr' or 'dense'.
    

    Returns
    --------------------
    C - hodlr | double
        The matrix of product. 


hlu.m
^^^^^^^^^^^^^

.. code:: matlab

        Compute LU factorization for HODLR matrix H.
    
        Parameters
        --------------------
        H - hodlr
            Matrix in HODLR format - hodlr class.
        
        oformat - str, default='hodlr'
            The output format. 'dense' or 'hodlr'.
    
        epsilon - double, default is the vareps of holdlr matrix H
            The vareps for recompression.
    
        Returns
        --------------------
        L - double
            The upper triangular matrix L is computed such that L * U = H. 
        U - double
            The upper triangular matrix U is computed such that L * U = H. 




mhlu.m
^^^^^^^^^^^^^

.. code:: matlab

    Compute LU factorization for HODLR matrix H.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    prec - precision
        Precision to simulate the factorization.

    oformat - str, default='hodlr'
        The output format. 'dense' or 'hodlr'.

    vareps - double, default is the vareps of holdlr matrix H
        The vareps for recompression.

    Returns
    --------------------
    L - double
        The upper triangular matrix L is computed such that L * U = H. 
    U - double
        The upper triangular matrix U is computed such that L * U = H. 






hchol.m
^^^^^^^^^^^^^

.. code:: matlab

    Compute Cholesky factorization for symmetric positive-definite HODLR matrix H.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    oformat - str, default='hodlr'
        The output format, either 'hodlr' or ''dense.


    Returns
    --------------------
    R - double
        The upper triangular matrix R is computed such that R' * R = H.  



mhchol.m
^^^^^^^^^^^^^


.. code:: matlab

    Compute Cholesky factorization for symmetric positive-definite HODLR matrix H.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    prec - precision
        Precision to simulate the factorization.

    oformat - str, default='hodlr'
        The output format, either 'hodlr' or ''dense.


    Returns
    --------------------
    R - double
        The upper triangular matrix R is computed such that R' * R = H.  



recover.m
^^^^^^^^^^^^^

.. code:: matlab

    The function is to recover a HODLR format into array format

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
 
    issparse - boolean
        `1` indicates returning sparse format, `0` indicates returning full arrary. 


    Returns
    --------------------
    A - double 
        Array in sparse or not.




hsize.m
^^^^^^^^^^^^^

.. code:: matlab

    The function is to return the shape of HOLDR matrix.

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    oformat - int, default=1
        If input is not the leafnode of the HODLR matrix:
            1: Outputs of varying number [m, n, m1, m2, n1, n2] for hierarchical block matrices.
            2: Outputs of varying number [m1, m2, n1, n2] for hierarchical block matrices.
        
        Otherwise:
            Outputs of varying number [m, n] for hierarchical block matrices.

    Returns
    --------------------
    [m, n, su1, su2, sv1, sv2] - int
        Indicates the size of rows and columns for H, size(H.U1, 1), size(H.U2, 1), size(H.V1, 2), size(H.V2, 2), respectively.




inverse.m
^^^^^^^^^^^^^

.. code:: matlab

    Parameters
    --------------------
    H - hodlr
        Matrix in HODLR format - hodlr class.
    
    algorithm - int, default=1
        The algorithm to implement inverse
    
    oformat - str, default = 'hodlr'
        The format of returns.
    
    Returns
    --------------------
    C - hodlr | double
        Return matrix in hodlr class or dense array.



lu_solve.m
^^^^^^^^^^^^^

.. code:: matlab

    Compute Hx = b using LU factorization.

    Parameters
    --------------------
    Setting 1:
        H - hodlr
            Matrix in HODLR format - hodlr class.

        b - double/single


    Setting 2:
        l - hodlr
            Matrix in HODLR format - hodlr class.

        U - hodlr
            Matrix in HODLR format - hodlr class.
            
        b - double/single

    Returns
    --------------------
    x - double
        The solution.





hstorage.m
^^^^^^^^^^^^^

.. code:: matlab

    Measure the theoretical storage of HODLR matrix

    Parameters
    --------------------
    
    H - hodlr, mphodlr, and amphodlr
        The input of HODLR matrix.


    Returns
    --------------------
    y - double 
        The number of bits for storage.