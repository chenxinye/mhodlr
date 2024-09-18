function [Q, R] = hqr(H, method)
    if strcmp(method, 'lintner')
        [Q, R] = lintner_qr(H);
    else % 'bebendorf')
        [Q, R] = bebendorf_qr(H);
    end
end

