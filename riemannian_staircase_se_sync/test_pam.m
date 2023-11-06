% A = [[1,2,3];[3,2,1];[2,1,3]];
A = [[4 0 1]; [2 3 2]; [1 0 4]];
% A = [[-11 -19 14]; [-6 -8 6]; [-12 -22 15]];


[eigvecs, eigvals] = eig(A)

thresh = 1e-5;

v_max = pim(A, thresh);

disp("v_max")
disp(v_max)

lambda_max = rayleigh_quotient(v_max, A); %Rayleigh quotient

disp("lambda_max")
disp(lambda_max)


A2 = A - (lambda_max) * eye(size(A));

v_max2 = pim(A2, thresh);

disp("v_max2")
disp(v_max2)

lambda_max2 = rayleigh_quotient(v_max2, A2); %Rayleigh quotient

disp("lambda_max2")
disp(lambda_max2)

disp("lambda_max2 + lambda_max")
disp(lambda_max2 + lambda_max)

