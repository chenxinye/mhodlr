classdef init_hodlr < hodlr
%{
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

%}
    properties(Access=public)

        U1 {mustBeNonNan, mustBeFinite, mustBeNumeric}
        V2 {mustBeNonNan, mustBeFinite, mustBeNumeric}

        U2 {mustBeNonNan, mustBeFinite, mustBeNumeric}
        V1 {mustBeNonNan, mustBeFinite, mustBeNumeric}
           
        D {mustBeNonNan, mustBeFinite, mustBeNumeric}

        A11 % HODLR format
        A22 % HODLR format
        
        level {mustBeInteger} = 0
        shape {mustBeNumeric} = []
        max_level {mustBeInteger} = 20
        bottom_level {mustBeInteger} = 0
        vareps {mustBeNonNan, mustBeFinite, mustBeNumeric} = 1.0e-12
        min_block_size {mustBeInteger} = 20
    end

    properties(Access=private)
        method {mustBeText} = 'svd'
        trun_norm_tp = '2'
    end

    methods(Access=public)
        function obj = init_hodlr(varargin)
            
        end
    end
end