%Compute the scalar algebraic error of the epipolar constraint and its derivatives
%function [e,JRT,HRT]=epipolarConstraint(R,T,x1,x2)
function [e,JRT,HRT]=epipolarConstraint(R,T,x1,x2)
Nx=size(x1,2);
x1=[x1;ones(1,Nx)];
x2=[x2;ones(1,Nx)];

flagComputeJacobian=false;
flagComputeSecondJacobian=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end

e=zeros(Nx,1);
E=epipolarBuildEFromRT(R,T);
for ix=1:Nx
    e(ix)=x1(:,ix)'*E*x2(:,ix);
end
if flagComputeJacobian
    JRT=zeros(1,6,Nx);
    if flagComputeSecondJacobian
        HRT=zeros(6,6,Nx);
    end
    for ix=1:Nx
        xL=x1(:,ix);
        xR=x2(:,ix);
        
        JRT(:,:,ix)=[-xL'*hat(T)*R*hat(xR) xR'*R'*hat(xL)];
        
        H1=-hat(R'*hat(T)*xL)*hat(xR);
        H2=hat(xR)*R'*hat(xL);
        HRT(:,:,ix)=[H1 H2;H2' zeros(3,3)];
    end
end