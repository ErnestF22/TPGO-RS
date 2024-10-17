%Compute the residual and derivatives of the rotation-only trifocal constraint
%function [e,ge,Dge]=trifocalRotationResiduals(n0,n1,n2,R1,R2)
%Compute the residual e=n0'*hat(R1*n1)*R2*n2, its gradient and its
%differential of the gradient w.r.t. R1,R2.
%Inputs
%   n0,n1,n2    [3xN] arrays of N normals in \real{3}
%   R1,R2       [3x3] rotations
%Outputs
%   e           [1xN] array where e(i)=n0(:,i)'*hat(R1*n1(:,i))*R2*n2(:,i)
%   ge          [6xN] array where ge(:,i) is the gradient of e(i)
%   Dge         [6x6xN] array where Dge(:,:,i) is the differential of ge(:,i)
function [e,ge,Dge]=trifocalRotationResiduals(R1,R2,n0,n1,n2)
Nn=size(n0,2);
flagComputeJacobian=false;
flagComputeSecondJacobian=false;
if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end


e=zeros(1,Nn);
if flagComputeJacobian
    ge=zeros(6,Nn);
    if flagComputeSecondJacobian
        Dge=zeros(6,6,Nn);
    end
end

for in=1:Nn
    n0i=n0(:,in);
    n1i=n1(:,in);
    n2i=n2(:,in);
    
    e(in)=n0i'*hat(R1*n1i)*R2*n2i;
    if flagComputeJacobian
        ge(:,in)=[hat(n1i)*R1'*hat(R2*n2i)*n0i; -hat(n2i)*R2'*hat(R1*n1i)*n0i];
        if flagComputeSecondJacobian
            Dge(:,:,in)=[
                 hat(n1i)*hat(R1'*hat(R2*n2i)*n0i)  hat(n1i)*R1'*hat(n0i)*R2*hat(n2i);
                -hat(n2i)*R2'*hat(n0i)*R1*hat(n1i) -hat(n2i)*hat(R2'*hat(R1*n1i)*n0i)];
        end
    end
end
