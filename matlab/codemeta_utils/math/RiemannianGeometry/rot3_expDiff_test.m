function rot3_expDiff_test
resetRands()
methods={'closedForm','svd'};

%pick a rotation
R=rot_randn();
%pick a curve (line) x(t) in the tangent space at R, T_R SO(3)
[xVec,dxVec]=real_randGeodFun(randn(3,1));

x=@(t) rot_hat(R,xVec(t));
dx=@(t) rot_hat(R,dxVec(t));

Exp=@(x) rot_exp(R,x);
for iMethod=1:length(methods)
    DiffExp=@(x,dx) rot3_expDiff(R,x,dx,...
        'tangent','method',methods{iMethod});
    figure(iMethod) 
    funCheckDifferential(Exp,x,DiffExp,dx)
end
