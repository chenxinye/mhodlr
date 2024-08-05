API Reference
======================================

@precision 
-----------

precision.m

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
