%Compute the double differential of the normalization operation
%function ddxNormalized=cnormalizeDDiff(x,dx,ddx)
%Compute the d^2/dt^2 cnormalize(x) from x, d/dt x, and d^2/dt^2 x

%See math/notes/householder.tex for detailed derivations
function ddxNormalized=cnormalizeDDiff(x,dx,ddx)
d=size(x,1);
[xp,nx]=cnormalize(x);
[dxp,Dxp]=cnormalizeDiff(x,dx);
if nx==0
    ddxNormalized=zeros(d,1);
else
    ddxNormalized=-(2*dxp*xp'+xp*dxp')*dx/nx+Dxp*ddx;
end
