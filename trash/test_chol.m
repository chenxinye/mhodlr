function verbose = test_chol(varargin)
    if nargin == 0
        eps = 1.0e-8;
    else
        eps = varargin{1};
    end

    n = 10;
    A = spdiags(ones(n, 1) * [2 17 2],  -1:1, n, n);

    verbose = 1;
    hA = hodlr(A);
    R = hchol(hA, eps);

    if norm(full(R'*R - A), 2) > eps
        verbose = 0;
    end

end