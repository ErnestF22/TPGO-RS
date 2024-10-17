function homographyContinuousEstimateRefine_test
%resetRands(2)
NPlanes=2;
NFrames=3;

n=randn(3,NPlanes);
n(:,1)=cnormalize(n(:,1));
%make sure that planes are in front of the camera
n=n.*([1;1;1]*sign(n(3,:)));

v=randn(3,NFrames);
w=0.1*randn(3,NFrames);

NPairsTotal=NPlanes*NFrames;
H=zeros(3,3,NPairsTotal);
E=zeros(NPairsTotal,2);

cnt=1;
for iPlane=1:NPlanes
    for iFrame=1:NFrames
        H(:,:,cnt)=homographyContinuousFromWVN(w(:,iFrame),v(:,iFrame),n(:,iPlane));
        E(cnt,:)=[iFrame iPlane];
        cnt=cnt+1;
    end
end

sigma=0.3;
wInit=w+sigma*randn(size(w));
vInit=v+sigma*randn(size(v));
nInit=n+sigma*randn(size(n));

[wEst,vEst,nEst]=homographyContinuousEstimateRefine(H,E,wInit,vInit,nInit,'collect');
fprintf('Finished in %d iterations\n',size(wEst,3)-1)
c=homographyContinuousEstimateRefine_cost(H,E,wEst,vEst,nEst);
semilogy(c)

