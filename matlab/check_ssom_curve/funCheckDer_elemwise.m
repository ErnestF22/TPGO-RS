function funCheckDer_elemwise

A0 = rand(4,3);
dA0 = rand(4,3);

[A, dA, A0, dA0, ddA] = real_geodFun(A0, dA0);

disp('A0')
disp(A0)
disp('dA0')
disp(dA0)

disp("funCheckDer(A, dA)")
funCheckDer(A, dA)

disp("funCheckDer(dA, ddA)")
funCheckDer(dA, ddA)

B = rand(5,4);
C = rand(3, 6);

A2 = @(t)B * A(t) * C;
dA2 = @(t)B * dA(t) * C;

disp("funCheckDer(B*A*C, B*dA*C)")
funCheckDer(A2, dA2)

end %file function