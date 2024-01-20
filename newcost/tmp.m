d = 3;
N = 2;


% A = randi([-10, 10],3, 3,5,5);
A = randi([-10, 10], d*N, d*N);
disp("A");
disp(A);

A_1 = vec(A');
disp("A_1");
disp(A_1);

A_2 = reshape(A_1, d*N*d, []);
disp("A_2");
disp(A_2);

A_3 = sum(A_2, 2);
disp("A_3");
disp(A_3);

A_4 = reshape(A_3, d*N,d);
disp("A_4");
disp(A_4);

A_out = A_4';
disp("A_out");
disp(A_out);

