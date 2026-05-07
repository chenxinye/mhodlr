========================================
Errata and Corrections
========================================

This page documents known issues and corrections in the mhodlr software package and related publications.


----

Issue 1: Variable Name Reuse in hstorage.m Creates Confusion
==============================================================

**Status:** Code refactored in version 3.1 (May 2026)

**Discovered by:** Ritesh Khan (Charles University)

**Date Reported:** November 2, 2025

**Severity:** Low - Code quality issue (does not affect correctness)

Description
-----------

A code quality issue was identified in ``hstorage.m`` where variable names are reused in a way that obscures the code's intent. Specifically, the variable ``k`` is used to store both the column count of ``U1`` and then overwritten with the row count of ``V2``.

**Affected Code (lines 44-47, 49-52, and similar in other sections):**

.. code-block:: matlab

    [m, k] = size(obj.U1);   % k = number of columns in U1
    [k, n] = size(obj.V2);   % k is OVERWRITTEN with number of rows in V2
    
    y = y + (m + n) * k * bits;  % Uses overwritten k

Technical Details
-----------------

**Affected Files:**

- ``mhodlr/hstorage.m`` (lines 44-52, 76-84, 108-116)

**Root Cause:**

The code reuses the variable ``k`` for two different purposes:

1. First, to store ``size(obj.U1, 2)`` - the number of columns in U1
2. Then, to store ``size(obj.V2, 1)`` - the number of rows in V2

For the factorization ``A12 = U1 * V2`` to be valid:

- If ``A12`` is of size ``m1 × n2``
- Then ``U1`` must be ``m1 × k`` 
- And ``V2`` must be ``k × n2``

The number of columns in ``U1`` **must equal** the number of rows in ``V2`` for matrix multiplication to be valid. Therefore, the overwriting of ``k`` does not affect numerical correctness, but it makes the code harder to understand and maintain.

Why This Does Not Affect Correctness
-------------------------------------

Due to the mathematical constraint of matrix multiplication, the storage calculation is correct:

.. code-block:: matlab

    % Current code computes:
    storage = (m + n) * k
    
    % Where m = rows of U1, n = cols of V2, k = rows of V2
    % This equals: m*k + n*k
    
    % Correct storage formula:
    storage = (rows_U1 * cols_U1) + (rows_V2 * cols_V2)
            = m * k1 + k2 * n
    
    % Since k1 == k2 (required for U1*V2 to be valid):
    storage = m * k + k * n = (m + n) * k  ✓

However, this equivalence is **implicit** and relies on mathematical properties that are not obvious from reading the code.

Impact
------

**Software Impact:**

- **No numerical errors**: Storage calculations are mathematically correct
- **Reduced code clarity**: Future maintainers may be confused by variable reuse
- **Maintenance risk**: Changes to the code could inadvertently break the implicit assumption


Suggested Improvement
---------------------

Ritesh Khan suggested a more explicit formulation:

.. code-block:: matlab

    [m, k1] = size(obj.U1);
    [k2, n] = size(obj.V2);
    
    y = y + (m*k1 + n*k2) * bits;

This makes the calculation explicit and would catch any unexpected cases where ``k1 ≠ k2`` (which would indicate a bug elsewhere).

Resolution
----------

**Code Improvement:**

Version 3.1 refactors ``hstorage.m`` to use more descriptive variable names:

.. code-block:: matlab

    % Improved version (v1.1.0)
    [m_U1, k_U1] = size(obj.U1);
    [k_V2, n_V2] = size(obj.V2);
    
    % Add assertion to verify matrix multiplication compatibility
    assert(k_U1 == k_V2, 'Inconsistent ranks: U1 and V2 cannot be multiplied');
    
    % Calculate storage explicitly
    storage_U1 = m_U1 * k_U1;
    storage_V2 = k_V2 * n_V2;
    y = y + (storage_U1 + storage_V2) * bits;


Acknowledgments
---------------

We thank **Ritesh Khan** from Charles University for:

- Carefully reviewing the code and identifying this subtle issue
- Suggesting a clearer formulation
- Contributing to improved code quality

This type of careful code review helps ensure the software is maintainable and understandable by the broader community.


Further Information
-------------------

For questions or concerns regarding these issues, please contact:

- **Email:** xinyechenai@gmail.com
- **GitHub Issues:** https://github.com/chenxinye/mhodlr/issues
- **Documentation:** https://mhodlr.readthedocs.io


Issue 2: Incorrect Storage Calculation for Adaptive Precision HODLR (amphodlr)
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
