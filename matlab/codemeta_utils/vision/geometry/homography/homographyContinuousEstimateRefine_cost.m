function c=homographyContinuousEstimateRefine_cost(H,E,w,v,n)
NIt=size(w,3);
c=zeros(NIt,1);
for it=1:NIt
    c(it)=evaluateCost(H,E,w(:,:,it),v(:,:,it),n(:,:,it));
end

function c=evaluateCost(H,E,w,v,n)
NPairs=size(H,3);

c=0;
for iPair=1:NPairs
    iFrame=E(iPair,1);
    iPlane=E(iPair,2);
    c=c+norm(H(:,:,iPair)-homographyContinuousFromWVN(w(:,iFrame),v(:,iFrame),n(:,iPlane)),'fro')^2;
end