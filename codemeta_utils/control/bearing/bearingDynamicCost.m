%Compute the total dynamic bearing cost and its derivative
%function [c,dc]=bearingDynamicCost(x,dx,yg,funs,alpha)
%Inputs
%   x       robot's current location
%   dx      derivative of x
%   yg      goal's bearing vectors
%   funs    reshaping function for the cost
%   alpha   weight of gradient component
%The derivative is computed with respect to the control given by 
%bearingDynamicControlDirect.
function [c,dc]=bearingDynamicCost(dx,y,ny,yg,funs,alpha)
[c1,gc1]=bearingCostGeneral(y,yg,ny,funs);
c2=(dx'*dx)/2;
c=alpha*c1+c2;

if nargout>1
    dY=bearingComputeDerivative(dx,y,ny);
    u=bearingDynamicControlDirect(y,dY,yg,funs,alpha);
    dc1=gc1'*dx;
    dc2=dx'*u;
    dc=alpha*dc1+dc2;
end
