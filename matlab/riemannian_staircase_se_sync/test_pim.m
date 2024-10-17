% A = [[1,2,3];[3,2,1];[2,1,3]];
A = [[4 0 1]; [2 3 2]; [1 0 4]];
% A = [[-11 -19 14]; [-6 -8 6]; [-12 -22 15]];


[eigvecs, eigvals] = eig(A)

thresh = 1e-5;

[lambda_max, v_max] = pim(A, thresh);

disp("v_max")
disp(v_max)

disp("lambda_max")
disp(lambda_max)


A2 = A - (lambda_max) * eye(size(A));

[lambda_max2, v_max2] = pim(A2, thresh);

disp("v_max2")
disp(v_max2)

disp("lambda_max2")
disp(lambda_max2)

disp("lambda_max2 + lambda_max")
disp(lambda_max2 + lambda_max)

