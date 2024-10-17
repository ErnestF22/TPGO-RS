function H=homographyContinuousEstimateLinear(x,dx)
NFrames=size(dx,3);
if NFrames>1
    H=zeros(3,3,NFrames);
    for iFrame=1:NFrames
        H(:,:,iFrame)=homographyContinuousEstimateLinear(x,dx(:,:,iFrame));
    end
else
    [A,b]=homographyContinuousEstimateLinearSystem(x,dx);
    [U,S,V]=svd(A,'econ');
    HL=reshape(V(:,1:8)*(diag(1./diag(S(1:8,1:8)))*(U(:,1:8)'*b)),3,3);
    l=svd(multisym(HL));
    H=HL-l(end)*eye(3);
end