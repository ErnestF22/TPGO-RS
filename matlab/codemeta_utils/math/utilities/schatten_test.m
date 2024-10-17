function schatten_test
p=0.7;
A=randn(5);

disp([schatten(A,p)^p trace((A*A')^(p/2))])
