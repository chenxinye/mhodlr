.. CLASSIX documentation

Welcome to MHODLR's documentation!
===================================


HODLR matrices are formulated by hierarchically partitioning the matrix in terms of a binary cluster tree and all off-diagonal blocks of each level of the tree are represented as low-rank matrices. This repository is concerned with Hierarchical Off-Diagonal Low-Rank (HODLR) matrices; we implement HODLR computations in MATLAB, and aim to provide a convenient API for HODLR operations. We also provide mixed precision simulation code for HODLR matrix computing.   


.. image:: demo.png
    :width: 640


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
