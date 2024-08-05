Build HODLR matrix
======================================

Now we try to build HODLR matrix of depth 5 and minmum block size of 2, using the SVD (singular value decomposition) truncation with approximation error of 1e-4; the normal HODLR format can be built via class ``@hodlr`` via:

.. code:: matlab

   rng(0); %fix randomness
   A = rand(500);
   depth = 5;
   min_block_size = 2;
   hA = hodlr(A, depth, min_block_size, 'svd', epsilon); % or simply use ``hA = hodlr(A)`` by omitting other parameters as default

The generators U, V (we will detail later) are represented in sparse format, the level indicates the current level of the HODLR matrix. 


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

  hA = 

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




The output refers to the attributes of HODLR matrix. The U, V represents generators that represent the off-diagonal blocks, i.e., :math:`A_{12} \approx U_1 V_2` and :math:`A_{21} \approx U_2 V_1`. 
The D is referred to as the diagonal block if the layer, indicated by the ``level``, is the bottom layer.  ``max_level`` refers to the maximum depth the HODLR matrix should have, but the practical depth is indicated in ``bottom_level``.
The ``vareps`` is the truncation error defined by the users (as input). The ``min_block_size`` is set to 2 (as input) as mentioned above . 

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