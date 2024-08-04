
n = 1000;
A = rand(n, n); % generate test matrix
hA = hodlr(A, 5, 100);
[L, U] = hlu(hA);

norm(hdot(L, U, 'dense') - A, 'fro')


n = 100;

A = rand(n, n);
disp(cond(A, 2));

%$ working precision
epsilon = 0.00001;
hA = hodlr(A, 5, 20, 'svd', epsilon);

[L, U] = hlu(hA);

% hh = hdot(hA.transpose(), hA, 'hodlr')
% R = hchol(hh, 'hodlr');
% htrsu(hA, R, 1);
% rR = recover(R);
% norm(rR' * rR - recover(hh),2)
% RA = recover(hA);
% disp(norm(RA - A, 'fro'));
disp('Working precision, eps=0.001')
P = hpcond(hA, epsilon); % input a HODLR matrix
disp(cond(P*A, 2)); % condition number is about 1

disp('Working precision, eps=0.01')
epsilon = 0.0001;
P = hpcond(A, epsilon); % input a default matrix which will be automatically 
                        % converted into HODLR matrix
                        
disp(cond(P*A, 2)); % condition number decrease


%% mix precision
u1 = precision('q43');
u2 = precision('q43');
u3 = precision('q43');
u4 = precision('h');
u5 = precision('h');
u_chain = prec_chain(u1, u2, u3, u4, u5);

disp('Precision [q43, q43, q43, h, h], eps=0.001')
epsilon = 0.001;
hA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
% RA = recover(hA);
% disp(norm(RA - A, 'fro'));
P = mphpcond(hA, u_chain, epsilon); % input a HODLR matrix
disp(cond(P*A, 2)); % condition number is about 1

disp('Precision [q43, q43, q43, h, h], eps=0.01')
epsilon = 0.01;
hA = mphodlr(u_chain, A, 5, 20, 'svd', epsilon);
P = mphpcond(hA, u_chain, epsilon); % input a HODLR matrix
disp(cond(P*A, 2)); % condition number decrease

