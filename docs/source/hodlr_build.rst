Build HODLR matrix
======================================

Now we try to build HODLR matrix of depth 5 and minmum block size of 2, using the SVD (singular value decomposition) truncation with approximation error of 1e-4; the normal HODLR format can be built via class ``@hodlr`` via:

.. code:: matlab

   rng(0); %fix randomness
   A = rand(500);
   depth = 5;
   min_block_size = 2;
   hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

The generators U, V (U1, V1, U2, and V2; we will detail later) are represented in sparse format, the level indicates the current level of the HODLR matrix. 


The simulation to mixed precision HODLR matrix construction is implemented by the class ``@mphodlr`` and ``@amphodlr``.
Previous to the settings, users are required to provide the precisions (as specified in Section `Customized precision <https://mhodlr.readthedocs.io/en/stable/precision.html>`_.).
The following example illustrates the usage:

.. code:: matlab

   % define the precisions
   u1 = precision('d');
   u2 = precision('s');
   u3 = precision('h');
   u4 = precision('b');
   u5 = precision('q52');

   % build the collection of precisions
   u_chain = prec_chain(u1, u2, u3, u4, u5); % the order matters!

   % build the precisions according to u_chain (each layer uses the corresponding preicison)
   mphA = mphodlr(u_chain, A, depth, min_block_size, 'svd', epsilon); 
   mprA = recover(mphA); % recover from the HODLR format
   
   % build the precisions automatically
   aphA = amphodlr(u_chain, A, depth, min_block_size, 'svd', epsilon); 
   aprA = recover(aphA); % recover from the HODLR format



Now we plain the HODLR format. First we need print out the variables ``hA``, ``mphA``, and ``aphA``: 

.. code:: matlab

   disp(hA)


Result is 

.. code:: matlab

  hodlr with properties:

                U1: [250×249 double]
                V2: [249×250 double]
                U2: [250×249 double]
                V1: [249×250 double]
                 D: []
               A11: [1×1 hodlr]
               A22: [1×1 hodlr]
             level: 1
             shape: [500 500]
         max_level: 5
      bottom_level: 5
            vareps: 1.0000e-04
    min_block_size: 2




The output refers to the attributes of HODLR matrix. The U, V (U1, V1, U2, and V2) represent generators that represent the off-diagonal blocks, i.e., :math:`A_{12} \approx U_1 V_2` and :math:`A_{21} \approx U_2 V_1`. 
The A11 and 
The D is referred to as the diagonal block if the layer, indicated by the ``level``, is the bottom layer.  ``max_level`` refers to the maximum depth the HODLR matrix should have, but the practical depth is indicated in ``bottom_level``.
The ``vareps`` is the truncation error defined by the users (as input). The ``min_block_size`` is set to 2 (as input) as mentioned above . 
The ``shape`` refers to the size of the dense matrix.  

.. code:: matlab

   disp(mphA)


Result is 

.. code:: matlab

  mphodlr with properties:

                U1: [250×249 double]
                V2: [249×250 double]
                U2: [250×249 double]
                V1: [249×250 double]
                 D: []
               A11: [1×1 mphodlr]
               A22: [1×1 mphodlr]
             level: 1
     prec_settings: {1×5 cell}
             shape: [500 500]
         max_level: 5
      bottom_level: 5
            vareps: 1.0000e-04
    min_block_size: 2


The mhodlr object contains an additional parameter ``prec_settings``, which indicates the precision used in each layer. 

.. admonition:: Note

   If the size of the collection of precisions is less than the depth, the rest of the layers will use double precision, as indicated in the warning information.

   .. code:: 

      Warning: The number of precisions used are less than the maximum
      tree level that can achieve. The remaining level will use the
      working precision for compresion. 


Similarly, by printing out the ``aphA``, we get 

.. code:: matlab

  amphodlr with properties:

                U1: [250×249 double]
                V2: [249×250 double]
                U2: [250×249 double]
                V1: [249×250 double]
                 D: []
               A11: [1×1 amphodlr]
               A22: [1×1 amphodlr]
             level: 1
             shape: [500 500]
         max_level: 5
      bottom_level: 5
         normOrder: [8.3423e+04 2.0891e+04 … ] (1×6 double)
         precIndex: [2 2 2 2 2]
      unitRoundOff: [1.1102e-16 1.1102e-16 … ] (1×6 double)
    min_block_size: 2
            vareps: 1.0000e-04
     prec_settings: {1×6 cell}


As shown in the output, we get three more parameters\: ``normOrder``, ``precIndex`` and ``unitRoundOff`` which separately denote the norm value of each layer, the precision used in each layer (indicated by the order of u_chain) and the unit-roundoff of each precision, respectively. Note we got six elements for  ``normOrder`` and ``unitRoundOff``, which more than the depths since the first element of them correspond the the layer 0. 