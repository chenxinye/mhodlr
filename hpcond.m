function P = hpcond(A, epsilon, max_level, min_block_size, method)
    
    if nargin < 2
        epsilon = 1e-5;
        max_level = 5;
        min_block_size = 50;
        method = 'svd';
    elseif nargin == 2
        max_level = 5;
        min_block_size = 50;
        method = 'svd';
    end

    if ~isa(A, 'hodlr')
        N = size(A, 1);
        A = hodlr(A, max_level, min_block_size, method, epsilon);
    else
        N = hsize(A, 1);
    end
    
    P = lu_solve(A, eye(N));
end