%function lambda=homographyTriangulateDepths(H,x1,x2)
%Compute lambdas such that lambda(2)*[x2;1]=lambda(1)*H*[x1;1], lambda(1)=1 
function lambda=homographyTriangulateDepths(H,x1,x2)
x1=homogeneous(x1);
x2=homogeneous(x2);

NPoints=size(x1,2);
lambda=zeros(2,NPoints);
for iPoint=1:NPoints
    [U,S,V]=svd([H*x1(:,iPoint) x2(:,iPoint)]);
    lambda(:,iPoint)=-V(:,2)/V(1,2);       %make lambda(1)=1
end
