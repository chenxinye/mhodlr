function [Q, R] = hqr(H, method)
    if strcmp(method, 'lintner')
        [Q, R] = lintner_qr(H);

    else if strcmp(method, 'bebendorf')
        [Q, R] = bebendorf_qr(H);

    else if strcmp(method, 'kressner')
        [Q, R] = kressner_qr(H);
    end
end

