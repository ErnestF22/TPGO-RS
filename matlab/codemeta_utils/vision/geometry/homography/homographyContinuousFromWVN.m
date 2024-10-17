function H=homographyContinuousFromWVN(w,v,nVec)
if size(nVec,1)==4
    nVec=planeNVecToNScaled(nVec);
end
N=size(w,2);
H=zeros(3,3,N);
for iN=1:N
    H(:,:,iN)=-hat3(w(:,iN))+v(:,iN)*nVec(:,iN)';
end