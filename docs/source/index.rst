.. MHODLR documentation

Welcome to MHODLR's documentation!
===================================

HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. This repository is concerned with Hierarchical Off-Diagonal Low-Rank (HODLR) matrices; we implement HODLR computations in Matlab and aim to provide a convenient API for HODLR operations. We also provide mixed precision simulation code for HODLR matrix computing.   

It is widely known that low precision can reduce data communication and be more energy- and storage-efficient. Regarding the IEEE standard for floating point, single precision arithmetic can be twice as fast as double precision on specific hardware and the half-precision arithmetic achieves 4 times speedup over double precision. Using the software ``mhodlr``, one can know what precisions are required for the HODLR matrix construction by evaluating their reconstruction error and computations error by simulating various precisions.  

.. image:: demo.png
    :width: 960


Guide
-------------

.. toctree::
   :maxdepth: 2
   
   start.rst
   hodlr_build.rst
   matrix_compute.rst
   precision.rst
   api_ref.rst

Others
-------------
.. toctree::
   :maxdepth: 2
   
   teams.rst
   acknow.rst
   license.rst


Indices and Tables
-------------------

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
