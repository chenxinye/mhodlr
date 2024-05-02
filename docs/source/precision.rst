Customized precision
======================================

With the software ``mhodlr``, you can easily customized your preferred precision for HODLR simulation. 
The precision is defined by class ``precision``, you can simply use a few letters to define the precision, or specify the detail of the precision used.

.. code:: matlab

  precision(6, 10, 1, 1, 0.5, 1)

.. code:: bash

  >> ans = 
  
    precision with properties:
  
              t: 6
           emax: 15
          round: 10
      subnormal: 1
         explim: 1
           prob: 0.5000
           flip: 1
       randfunc: @(n)rand(n,1)
