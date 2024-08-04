Customized precision
======================================

With the software ``mhodlr``, you can easily customized your preferred precisions for the HODLR format conversion. 
The precision is defined by class ``@precision``, you can simply use a few letters to define the precision, or specify the detail of the precision used.

We give a few examples to illustrate its usage. 

One can simply specify the common floating point format by entering a given name or abbreviation (for detail, see the API reference as below)

.. code:: matlab

  u = precision('s');


.. code:: bash

  >> ans = 
  precision with properties:
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

Then, the result is:

.. code:: bash

  >> ans = 
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


Or you can leave the ``emax`` empty to use default value for it, then it can be: 

.. code:: matlab

  precision(6, 2, 1, 1, 0.5, 1)

the result is similar:

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
