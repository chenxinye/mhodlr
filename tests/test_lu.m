function b = test_lu(varargin)
    if nargin == 0
        eps = 1.0e-8;
    else
        eps = varargin{1};
    end

    n = 10;
    A = spdiags(ones(n, 1) * [2 17 -11],  -1:1, n, n);

    b = 1;
    hA = hodlr(A);
    [L, U] = hlu(hA, eps);

    if norm(full(L*U - A), 2) > eps
        b = 0;
    end

end