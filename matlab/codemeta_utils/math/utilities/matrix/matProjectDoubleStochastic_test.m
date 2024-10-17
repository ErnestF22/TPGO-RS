function projectDoubleStochastic_test
%resetRands()
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
cvx_end

XB=projectDoubleStochastic(Z);

disp([X-XB])
