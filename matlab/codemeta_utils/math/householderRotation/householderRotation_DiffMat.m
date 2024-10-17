%Compute the matrix representation of the differential for householderRotations
%function DH=householderRotation_DiffMat(x1,x2)
function DH=householderRotation_DiffMat(x1,x2)
D1=size(x1,1);
D2=size(x2,1);
[N,x1,x2]=argRepmat(2,x1,x2);
flagStandardVector=(D2==1);

if N>1
    if flagStandardVector
        DH=zeros(D1,D1,N);
    else
        DH=zeros(D1,2*D1,N);
    end
    for iN=1:N
        DH(:,:,iN)=householderRotation_Diff(x1(:,iN),x2(:,iN));
    end
else
    [Dx1p,x1p]=cnormalizeDiffMat(x1);

    if flagStandardVector
        I=eye(D1);
        x2p=I(:,x2);
        Dx2p=[];
    else
        [Dx2p,x2p]=cnormalizeDiffMat(x2);
    end
        
    v=x1p+x2p;

    [Dvp,vp]=cnormalizeDiffMat(v);

    DH=-2*hat3(vp)*Dvp*[Dx1p Dx2p];
end
