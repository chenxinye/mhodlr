function verbose = test_dot(varargin)
    if nargin == 0
        eps = 1.0e-8;
    else
        eps = varargin{1};
    end

    M = rand(5,3,2);
    C = M(:,:,1) * M(:,:,2)';

    verbose = 1;

    hA = hodlr(M(:,:,1));
    hB = hodlr(M(:,:,2)');

    if norm(hdot(hA, M(:,:,2)', 'double') - C, 2) > eps
        verbose = 0;
    end

    if norm(hdot(hA, hB, 'double') - C, 2) > eps
        verbose = 0;
    end

    if norm(hdot(M(:,:,1), hB, 'double') - C, 2) > eps
        verbose = 0;
    end
end