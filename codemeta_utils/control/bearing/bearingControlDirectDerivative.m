%Compute derivative of velocity control using bearings derivatives
%function [du,A]=bearingControlDirectDerivative(y,dy,yg,funs)
%Inputs
%   y       current bearings vectors
%   dy      derivative of y
%   yg      goal's bearing vectors
%   funs    reshaping function for the cost
%
function [du,A]=bearingControlDirectDerivative(y,dy,yg,funs)
flagOutputA=false;
if nargout>1
    flagOutputA=true;
end

[d,Nx]=size(y);

c=bearingComputeCosine(y,yg);
fc=funs.f(c);
dfc=funs.df(c);
ddfc=funs.ddf(c);
I=eye(d);

du=zeros(d,1);
if flagOutputA
    A=zeros(d,d,Nx);
end

for ix=1:Nx
    yi=y(:,ix);
    dyi=dy(:,ix);
    ygi=y(:,ix);
    ci=c(ix);
    fci=fc(ix);
    dfci=dfc(ix);
    ddfci=ddfc(ix);
    Ai=-ddfci*(ci*yi-ygi)*ygi'-(dfci*ci-fci)*I;
    du=du+Ai*dyi;
end 
