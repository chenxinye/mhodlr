function nrm = hnorm(hA, varargin)
    if nargin == 1
        p = "fro";
    else
        p = varargin{1};
    end

    hA = recover(hA);
    nrm = norm(hA, p);
end