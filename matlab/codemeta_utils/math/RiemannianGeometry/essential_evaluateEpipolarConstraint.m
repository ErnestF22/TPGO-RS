function [e,JQ,HQ]=essential_evaluateEpipolarConstraint(Q,x1,x2)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end

Nx=size(x1,2);
x1=[x1;ones(1,Nx)];
x2=[x2;ones(1,Nx)];

E=essential_getE(Q);

e=sum(x1.*(E*x2));

if flagComputeJacobian
    JQ=zeros(6,Nx);
    if flagComputeSecondJacobian
        HQ=zeros(6,6,Nx);
    end
    for ix=1:Nx
        x1i=x1(:,ix);
        x2i=x2(:,ix);
        Etx1i=E'*x1i;
        Ex2i=E*x2i;
        hx1i=hat(x1i);
        hx2i=hat(x2i);
        
        
        JQ(1:3,ix)=hx1i*Ex2i;
        JQ(4:6,ix)=hx2i*Etx1i;
        if flagComputeSecondJacobian
            HQ(1:3,1:3,ix)=hx1i*hat(Ex2i);
            HQ(4:6,1:3,ix)=-hx2i*E'*hx1i;
            HQ(1:3,4:6,ix)=HQ(4:6,1:3,ix)';
            HQ(4:6,4:6,ix)=hx2i*hat(Etx1i);
        end
    end
end
