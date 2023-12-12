function metric_decomp_test
M=randn(2);
M=M'*M;
beta=rand;
alpha=randn(2,1);

MM=[beta alpha'; alpha M];

[U,betap]=metric_decomp(M,alpha,beta);
Uinv=inv(U);
disp(U'*MM*U-blkdiag(betap,M))
disp(MM-Uinv'*blkdiag(betap,M)*Uinv)
