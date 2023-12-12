function yCenter=bearingFitFov(y)
if size(y,2)==1
    y=permute(y,[1 3 2]);
end

Ny=size(y,2);
a=sphere2_distSigned(eye(2,1),y);
[aSorted,idxSort]=sort(a);

aSortedPair=aSorted-circshift(aSorted,[0,1]);
aSortedPair(1)=aSortedPair(1)+2*pi;


[aPairMax,idxPair]=max(aSortedPair);
idxyPair=[idxSort(idxPair) idxSort(mod(idxPair-2,Ny)+1)];
yCenter=cnormalize(sum(y(:,idxyPair),2));

if aPairMax<pi
    yCenter=-yCenter;
end

