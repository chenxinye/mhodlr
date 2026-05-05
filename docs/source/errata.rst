========================================
Errata and Corrections
========================================

This page documents known issues and corrections in the mhodlr software package and related publications.


----

Issue #1: Incorrect Storage Calculation for Adaptive Precision HODLR (amphodlr)
================================================================================

**Status:** Fixed in version 3.1 (May 2026)

**Date Reported:** April 28, 2026

**Severity:** High - Affects storage calculation accuracy

Description
-----------

A bug was identified in the ``hstorage.m`` function when computing storage requirements for adaptive precision HODLR matrices (``amphodlr`` class). The issue involves two related problems:

**Problem 1: Precision Index Propagation Failure**

The ``precIndex`` array, which stores the precision level used at each hierarchical level, was correctly computed for the root node during matrix construction but was not propagated to child nodes (``A11``, ``A22``, etc.). This caused ``hstorage`` to use incorrect (default) precision values when recursively calculating storage for the tree structure.

**Problem 2: Array Indexing Inconsistency**

There was an index offset mismatch between how ``precIndex`` values were stored in ``build_hodlr_mat`` (with ``-1`` offset) and how they were accessed in ``hstorage.m`` (without compensating ``+1``). This resulted in reading incorrect precision settings from the ``prec_settings`` array. So the storage result should be lower than previously presented.

Technical Details
-----------------

**Affected Files:**

- ``mhodlr/@amphodlr/amphodlr.m`` (lines 191-199, 349, 361-365)
- ``mhodlr/hstorage.m`` (line 103)

**Root Cause:**

In ``amphodlr.m`` constructor:

.. code-block:: matlab

    obj.precIndex = ones(1, obj.max_level);  % initialized
    [obj, obj.precIndex, obj.precIndexBool] = build_hodlr_mat(obj, A, obj.level, ...
                                                    obj.precIndex, obj.precIndexBool);
    obj.precIndex = obj.precIndex(1: obj.bottom_level);  % truncated to actual depth

The ``precIndex`` array was updated only for the root node. When ``build_hodlr_mat`` recursively created child nodes:

.. code-block:: matlab

    [obj.A11, precIndex, precIndexBool] = build_hodlr_mat(obj, A(1:rowSplit, 1:colSplit), ...
        level, precIndex, precIndexBool);

The returned ``precIndex`` array was not assigned back to ``obj.A11.precIndex``, leaving child nodes with the default initialization value of ``ones(1, max_level)``.

Additionally, in ``hstorage.m`` line 103:

.. code-block:: matlab

    bits = obj.prec_settings{obj.precIndex(level)}.bits;  % Missing +1!

Should have been:

.. code-block:: matlab

    bits = obj.prec_settings{obj.precIndex(level)+1}.bits;  % Correct indexing

to match the indexing convention used in ``build_hodlr_mat``.

Impact
------

**Software Impact:**

- The ``hstorage`` function **overestimated** storage requirements for ``amphodlr`` matrices
- Child nodes were assumed to use higher precision than actually used
- The error magnitude depends on tree depth and precision differences between levels

**Publication Impact:**

This bug affects the storage calculation results reported in:

    Chen, X., Carson, E., & Liu, X. (2024). Mixed precision HODLR matrices. 
    *SIAM Journal on Scientific Computing*, 46(3), A1408-A1435. 
    https://doi.org/10.1137/23M1546592

Verification
------------

Users can verify the fix by running the test suite:

.. code-block:: matlab

    addpath('mhodlr/tests')
    run_storage_tests  % New test suite for storage calculations

Or manually verify on a small example:

.. code-block:: matlab

    rng(0);
    A = rand(64, 64);
    u_chain = prec_chain(precision('d'), precision('s'), precision('h'));
    
    % Build amphodlr matrix
    aphA = amphodlr(u_chain, A, 5, 4, 'svd', 1e-8);
    
    % Calculate storage (now corrected)
    storage_bits = hstorage(aphA);
    
    % Verify child nodes have correct precIndex
    assert(~all(aphA.A11.precIndex == 1), 'precIndex should be propagated');

Acknowledgments
---------------

We sincerely thank **Jindřich Pohl** for:

- Discovering and thoroughly documenting this issue
- Providing a clean and effective fix
- Contributing to the quality and reliability of this software

His careful code review and contribution exemplify the best practices of open-source scientific computing.

Further Information
-------------------

For questions or concerns regarding this erratum, please contact:

- **Email:** xinyechenai@gmail.com
- **GitHub Issues:** https://github.com/chenxinye/mhodlr/issues
- **Documentation:** https://mhodlr.readthedocs.io

----

**Citation for Erratum:**

If you use this software in your research, please cite both the original paper and acknowledge this correction:

.. code-block:: bibtex

    @article{chen2024mixed,
      title={Mixed precision HODLR matrices},
      author={Chen, Xinye and Carson, Erin and Liu, Xiaobo},
      journal={SIAM Journal on Scientific Computing},
      volume={46},
      number={3},
      pages={A1408--A1435},
      year={2024},
      note={Erratum: Storage calculation bug fixed in software v1.1.0 (May 2026). 
            Corrected results show improved storage efficiency.}
    }

----

*This erratum reflects our commitment to transparency and scientific integrity. 
We appreciate the community's vigilance in ensuring the quality of scientific software.*
