%Return 3x3 Householder rotation which is closest to the identity
%function HMin=householderRotation3Min(x1,x2)
function HMin=householderRotation3Min(x1,x2)
H=householderRotation(x1,x2);

if size(x2,1)==1
    Rx2pi=diag(-ones(3,1));
    Rx2pi(x2,x2)=1;
else
    Rx2pi=2*(x2*x2')-eye(3);
end
HMin=H*Rx2pi;
