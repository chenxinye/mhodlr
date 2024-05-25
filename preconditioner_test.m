n = 100;

A = rand(n, n);
disp(cond(A, 2));

epsilon = 0.00001;
hA = hodlr(A, 5, 20, 'svd', epsilon);
% RA = recover(hA);
% disp(norm(RA - A, 'fro'));
P = hpcond(hA, epsilon); % input a HODLR matrix
disp(cond(P*A, 2)); % condition number is about 1

epsilon = 0.01;
P = hpcond(A, epsilon); % input a default matrix which will be automatically 
                        % converted into HODLR matrix
                        
disp(cond(P*A, 2)); % condition number decrease

