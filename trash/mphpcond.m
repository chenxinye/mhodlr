function P = mphpcond(A, u_chain, epsilon, max_level, min_block_size, method)
    
    if nargin < 3
        epsilon = 1e-5;
        max_level = 5;
        min_block_size = 20;
        method = 'svd';
    elseif nargin == 3
        max_level = 5;
        min_block_size = 20;
        method = 'svd';
    end

    if ~isa(A, 'mphodlr')
        N = size(A, 1);
        A = mphodlr(u_chain, A, max_level, min_block_size, method, epsilon);
    else
        N = hsize(A, 1);
    end
    
    P = lu_solve(A, eye(N));
end