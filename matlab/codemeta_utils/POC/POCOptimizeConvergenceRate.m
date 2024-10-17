function POCOptimizeConvergenceRate
resetRands()

S1=randn(3);
S2=randn(3);
A1=S1*diag([-1,-2,-3])/S1;
A2=S2*diag([1,-2,-3])/S2;
A3=S2*diag([-1,2.5,2.8])/S2;

epsilon=2.22;
cvx_begin sdp quiet
    %variable P(3,3) symmetric
    variables c2 c3
    c2>=0;
    c3>=0;
    A=A1+c2*A2+c3*A3;
    A'+A<=-epsilon*eye(3);
cvx_end
cvx_status
cvx_optval

disp(c2)
disp(c3)

disp(A)
disp(eig(A))
disp(eig(A+A'))
