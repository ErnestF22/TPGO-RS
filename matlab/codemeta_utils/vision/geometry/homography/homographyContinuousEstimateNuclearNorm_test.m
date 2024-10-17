function homographyContinuousEstimateNuclearNorm_test
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

%H=H+0.05*randn(size(H));

[wEst,vEst,nEst,output]=homographyContinuousEstimateNuclearNorm(H,E);

%compute ground truth
MOuterTruth=v(:)*n(:)';

disp('[MOuterTruth MOuter]')
disp([MOuterTruth output.MOuter])

disp('[H MOuter]')
compareHAndM(H,output.MOuter,E)

disp('[H MOutherTruth]')
compareHAndM(H,MOuterTruth,E)

disp('svd(MOuter)')
disp(svd(output.MOuter)')

disp('[v vEst]')
disp([v vEst])
disp('[n nEst]')
disp([n nEst])
disp('[w wEst]')
disp([w wEst])

end

function compareHAndM(H,M,E)
NFrames=max(E(:,1));
NPlanes=max(E(:,2));

idxFrames=reshape(1:3*NFrames,3,NFrames);
idxPlanes=reshape(1:3*NPlanes,3,NPlanes);

NPairs=size(H,3);

for iPair=1:NPairs
    iFrame=E(iPair,1);
    iPlane=E(iPair,2);
    disp([multisym(H(:,:,iPair)) multisym(M(idxFrames(:,iFrame),idxPlanes(:,iPlane)))])
end
end