function [A,b]=homFlowParametersEstimate4ptTranslationSystemMat(x,w)
NX=size(x,2);
A=zeros(2*NX,9);
b=zeros(2*NX,1);
idxA=reshape(1:2*NX,2,NX);
for iX=1:NX
    xi=x(1,iX);
    yi=x(2,iX);
    A(idxA(:,iX),:)=kron([eye(2) -x(:,iX)],[x(:,iX)' 1]);
        
    b(idxA(:,iX),:)=[-w(2)+yi*w(3)+xi*yi*w(1)-xi^2*w(2);
        w(1)-xi*w(3)+yi^2*w(1)-xi*yi*w(2)];
end
