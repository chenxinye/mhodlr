function verbose = test_recover(varargin)
    if nargin == 0
        eps = 1.0e-8;
    else
        eps = varargin{1};
    end

    A = rand(5,6,3);
    hA1 = hodlr(A(:, :, 1));
    hA2 = hodlr(A(:, :, 2), 'qr', 1.0e-8, 2);
    hA3 = hodlr(A(:, :, 3), 'qr', 1.0e-8, 2, 3);

    verbose = 1;
    if norm(recover(hA1) - A(:, :, 1), 2) > eps
        verbose = 0;
    end

    if norm(recover(hA2) - A(:, :, 2), 2) > eps
        verbose = 0;
    end

    if norm(recover(hA3) - A(:, :, 3), 2) > eps
        verbose = 0;
    end
end
