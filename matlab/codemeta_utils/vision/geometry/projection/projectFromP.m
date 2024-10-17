%Project points given projection matrices and compute Jacobians
%function [x,Jx]=projectFromP(P,X)
%
function [x,Jx,Hx]=projectFromP(P,X)
flagComputeJacobian=false;
flagComputeSecondJacobian=false;

if nargout>1
    flagComputeJacobian=true;
end
if nargout>2
    flagComputeSecondJacobian=true;
end

%get dimensions
NX=size(X,2);
NP=size(P,3);
d=size(P,1)-1;  %image space dimension
D=size(P,2)-1;  %domain space dimension

%initialize matrices
x=zeros(d,NX,NP);
if flagComputeJacobian
    Jx=zeros(d,D,NX,NP);
end
if flagComputeSecondJacobian
    Hx=zeros(D,D,d,NX,NP);
end

%for each view (i.e., projection matrix)
for iP=1:NP
    P1=P(:,:,iP);
    %intermediate product P*X
    PX=P1(:,1:end-1)*X+P1(:,end)*ones(1,NX);
    x(:,:,iP)=PX(1:2,:)./([1;1]*PX(3,:));
    
    if flagComputeJacobian
        %for each point and the current view
        for iX=1:NX
            Jx(:,:,iX,iP)=(PX(end,iX)*P1(1:end-1,1:end-1)-PX(1:end-1,iX)*P1(end,1:end-1))/(PX(end,iX)^2);
        end
    end
    if flagComputeSecondJacobian
        P1d1=P1(end,1:end-1);
        ed1PX=PX(end,iX);
        for iX=1:NX
            for k=1:d
                P1k=P1(k,1:end-1);
                ekPX=PX(k,iX);
                
                Hx(:,:,k,iX,iP)=(ed1PX*(P1d1'*P1k-P1k'*P1d1)-2*P1d1'*(ed1PX*P1k-ekPX*P1d1))/ed1PX^3;
            end
        end
    end
end
