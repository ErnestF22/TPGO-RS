function projectPSD_test
n=6;
Z=randn(n);

cvx_begin
    variable X(n,n)
    minimize norm(X-Z,'fro')
    subject to
        X==semidefinite(n);
cvx_end

XB=projectPSD(Z);

disp([X XB])
