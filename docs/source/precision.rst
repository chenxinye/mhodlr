Customized precision
======================================

With the software ``mhodlr``, you can easily customized your preferred precision for HODLR simulation. 
The precision is defined by class ``precision``, you can simply use a few letters to define the precision, or specify the detail of the precision used.

We give a few examples to illustrate its usage. 

One can simply specify the common floating point format by entering a given name or abbreviation (for detail, see the API reference as below)

.. code:: matlab

  u = precision('s');


.. code:: bash

  >> ans = 

  precision with properties:
  >> ans = 
            t: 24
         emax: 127
        round: 1
    subnormal: 1
       explim: 1
         prob: 0.5000
         flip: 0
     randfunc: @(n)rand(n,1)
            u: 5.9605e-08

To specify the tuple (t, emax), one can enter the first input as an array like
.. code:: matlab

  precision([6, 12], 2, 1, 1, 0.5, 1)

Then, the result is

.. code:: matlab

  precision with properties:

            t: 6
         emax: 12
        round: 2
    subnormal: 1
       explim: 1
         prob: 1
         flip: 0.5000
     randfunc: @(n)rand(n,1)
            u: 0.0156


Or you can leave the ``emax`` empty, use default value for it, then it can be 
.. code:: matlab

  precision(6, 2, 1, 1, 0.5, 1)


precision with properties:

.. code:: bash

  precision with properties:
  >> ans = 
            t: 6
         emax: 15
        round: 2
    subnormal: 1
       explim: 1
         prob: 1
         flip: 0.5000
     randfunc: @(n)rand(n,1)
            u: 0.0156


API Reference
---------------------------
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
