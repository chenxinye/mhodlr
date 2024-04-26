function C = blkoffdiag(varargin)
%% The function is to check whether HODLR matrix H is square or not.
%
% Parameters
% --------------------
% A, B, C ..., (variable number) - numeric
%      Inputs for block matrix
%
% Returns
% --------------------
% C - numeric
%     The off diagonal block matrices of C are [A, B, C, ..., (variable number)]
%   
    if nargin == 1 & length(size(varargin{1})) == 3
        [mi, ni, num] = size(varargin{1});
        ml = ones(num+1, 1);
        nl = ones(num+1, 1);

        for i = 1:num
            [mi, ni] = size(varargin{1}(:, :, i));
            ml(i+1)= mi; nl(i+1) = ni;
        end

        ml = cumsum(ml);
        nl = cumsum(nl);

        C = zeros(ml(end)-1, nl(end)-1);
        for i = 1:num
            C(ml(i): ml(i+1)-1, nl(end)-nl(i+1)+1:nl(end) - nl(i)) = varargin{1}(:, :, i);
        end

    else
        ml = ones(nargin+1, 1);
        nl = ones(nargin+1, 1);

        for i = 1:nargin
            [mi, ni] = size(varargin{i});
            ml(i+1)= mi; nl(i+1) = ni;
        end

        ml = cumsum(ml);
        nl = cumsum(nl);

        C = zeros(ml(end)-1, nl(end)-1);
        for i = 1:nargin
            C(ml(i): ml(i+1)-1, nl(end)-nl(i+1)+1:nl(end) - nl(i)) = varargin{i};
        end
    end
        
end