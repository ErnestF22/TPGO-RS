function multikron_test
A=randn(3,5,2);
B=randn(2,1);

disp(norm(vec(cat(3,kron(A(:,:,1),B),kron(A(:,:,2),B))-multikron(A,B)),inf))
