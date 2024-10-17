function R=bearingsToRotation2D(v1,v2)
N=size(v1,2);
R=zeros(2,2,N);
for iN=1:N
    cs=cnormalize([v2(1,iN) -v2(2,iN); v2(2,iN) v2(1,iN)]\v1(:,iN));
    c=cs(1);
    s=cs(2);
    R(:,:,iN)=[c -s; s c];
end

