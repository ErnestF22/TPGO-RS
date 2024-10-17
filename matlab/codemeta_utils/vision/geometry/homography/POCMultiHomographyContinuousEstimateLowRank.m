function POCMultiHomographyContinuousEstimateLowRank
resetRands()
NPlanes=2;
NFrames=3;
lambda=1;

n=randn(3,NPlanes);
n(:,1)=cnormalize(n(:,1));
v=randn(3,NFrames);
w=0.1*randn(3,NFrames);

NPairsTotal=NPlanes*NFrames;
H=zeros(3,3,NPairsTotal);
E=zeros(NPairsTotal,2);

cnt=1;
for iPlane=1:NPlanes
    for iFrame=1:NFrames
        H(:,:,cnt)=hat3(w(:,iFrame))+v(:,iFrame)*n(:,iPlane)';
        E(cnt,:)=[iFrame iPlane];
        cnt=cnt+1;
    end
end

NidxFrames=NFrames*3;
idxFrames=reshape(1:NidxFrames,3,NFrames);
idxPlanes=reshape((1:NPlanes*3)+NidxFrames,3,NPlanes);


NPairs=size(E,1);
K=3*(NPlanes+NFrames);

%method='penalty';
method='equality';

switch method
    case 'penalty'
        JPair=cell(NPairs,2);
        for iPair=1:NPairs
            iFrame=E(iPair,1);
            iPlane=E(iPair,2);
            J1=zeros(K,3);
            J2=zeros(K,3);
            J1(idxPlanes(:,iPlane),:)=eye(3);
            J2(idxFrames(:,iFrame),:)=eye(3);

            JPair{iPair,1}=J1;
            JPair{iPair,2}=J2;
        end

        cvx_begin
            variable MBig(K,K) symmetric semidefinite
            f=lambda*trace(MBig);
            for iPair=1:NPairs
                J1=JPair{iPair,1};
                J2=JPair{iPair,2};
                E=(H(:,:,iPair)+H(:,:,iPair)')-(J1'*MBig*J2+J2'*MBig*J1);
                f=f+E(:)'*E(:);
            end
            minimize f
        cvx_end
    case 'equality'
        cvx_begin
            variable MBig(K,K) symmetric semidefinite
            minimize trace(MBig)
            subject to
                for iPair=1:NPairs
                    iFrame=E(iPair,1);
                    iPlane=E(iPair,2);
                    H(:,:,iPair)+H(:,:,iPair)'==MBig(idxPlanes(:,iPlane),idxFrames(:,iFrame))+MBig(idxFrames(:,iFrame),idxPlanes(:,iPlane));
                end
                trace(MBig(idxPlanes(:,1),idxPlanes(:,1)))==1;
        cvx_end
end

%compute ground truth
mTruth=zeros(K,1);
for iFrame=1:NFrames
    mTruth(idxFrames(:,iFrame))=v(:,iFrame);
end
for iPlane=1:NPlanes
    mTruth(idxPlanes(:,iPlane))=n(:,iPlane);
end
MBigTruth=mTruth*mTruth';

disp('[H MBig]')
compareHAndM(H,MBig,idxPlanes,idxFrames)

disp('[H MBigTruth]')
compareHAndM(H,MBigTruth,idxPlanes,idxFrames)

disp('svd(MBig)')
disp(svd(MBig)')

disp('[mTruth m]')
[U,S,Ut]=svd(MBig,'econ');
m=U(:,1);
disp([mTruth m]')
keyboard

end

function compareHAndM(H,M,idxPlanes,idxFrames)
    cnt=1;
    for iPlane=1:size(idxPlanes,2)
        for iFrame=1:size(idxFrames,2)
            disp([multisym(H(:,:,cnt)) multisym(M(idxPlanes(:,iPlane),idxFrames(:,iFrame)))])
            cnt=cnt+1;
        end
    end
end