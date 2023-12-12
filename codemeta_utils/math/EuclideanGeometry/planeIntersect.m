%Compute the line at the intersection of pairs of planes
%function bOrthLine=planeIntersect(bOrthPlane1,bOrthPlane2)
function bOrthLine=planeIntersect(bOrthPlane1,bOrthPlane2)
NPlanes1=size(bOrthPlane1,2);
NPlanes2=size(bOrthPlane2,2);

bOrthLine=zeros(4,2,NPlanes1,NPlanes2);
for iPlane1=1:NPlanes1
    for iPlane2=1:NPlanes2
        [U,~,~]=svd([bOrthPlane1(:,iPlane1) bOrthPlane2(:,iPlane2)]);
        bOrthLine(:,:,iPlane1,iPlane2)=U(:,1:2);
    end
end
