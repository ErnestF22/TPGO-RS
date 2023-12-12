function matProjectPSDDoubleStochastic_test
n=6;
Z=randn(n);
%Z=(Z+Z')/2;

u=ones(n,1);
cvx_begin
    variable X(n,n)
    minimize norm(X-Z,'fro')
    subject to
        X*u==u;
        X'*u==u;
        X==semidefinite(n);
cvx_end

XB=matProjectPSDDoubleStochastic(Z);

disp([X-XB])

