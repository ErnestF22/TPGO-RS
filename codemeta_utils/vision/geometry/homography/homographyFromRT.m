function H=homographyFromRT(R,T,NVec)
[N,d]=planeNVecToNd(NVec);
NPoses=size(R,3);
NPlanes=size(NVec,2);
H=zeros(3,3,NPoses,NPlanes);
for iPose=1:NPoses
    for iPlane=1:NPlanes
        H(:,:,iPose,iPlane)=R(:,:,iPose)+T(:,iPose)*(N(:,iPlane)'/d(iPlane));
    end
end
H=squeeze(H);

