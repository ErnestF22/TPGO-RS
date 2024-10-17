%function [R,T,N,lambda]=homographyToRT(H,x1,x2)
%Decompose an homography matrix H into rotation and translation. The image
%points x1 and x2 are used to solve the twisted pair ambiguity
function [R,T,N,lambda,transfPairs]=homographyToRT(H,x1,x2)
if ~exist('x2','var') || isempty(x2)
    [x1,x2]=homographySplitData(x1);
end

NH=size(H,3);
if size(x2,3)~=NH
    error('size(x2,3) and size(H,3) must be the same.')
end

if NH>1
    R=zeros(3,3,NH);
    T=zeros(3,NH);
    N=zeros(3,NH);
    lambda=zeros(2,size(x1,2),NH);
    trasfPairs=cell(NH,1);
    for iN=1:NH
        [R(:,:,iN),T(:,iN),N(:,iN),lambda(:,:,iN),trasfPairs{iN}]=...
            homographyToRT(H(:,:,iN),x1,x2(:,:,iN));
    end
else
    H=homographyNormalize(H,x1,x2);

    [U,S,V] = svd(H'*H);

    if det(U) < 0
     U = -U;
    end

    s1 = S(1,1); s3 = S(3,3);
    v1 = U(:,1); v2 = U(:,2); v3 = U(:,3);
    u1 = (v1*sqrt(1-s3) + v3*sqrt(s1 -1))/sqrt(s1 - s3);
    u2 = (v1*sqrt(1-s3) - v3*sqrt(s1 -1))/sqrt(s1 - s3);

    U1 = [v2, u1, -hat3(v2)*u1];
    U2 = [v2, u2, -hat3(v2)*u2];
    W1 = [H*v2, H*u1, -hat3(H*v2)*H*u1];
    W2 = [H*v2, H*u2, -hat3(H*v2)*H*u2];

    N1 = -hat3(v2)*u1;
    N2 = -hat3(v2)*u2;

    transfPairs(1)=struct('R',W1*U1','T', (H - W1*U1')*N1,'N', N1);
    transfPairs(2)=struct('R',W2*U2','T', (H - W2*U2')*N2,'N', N2);
    transfPairs(3)=struct('R',W1*U1','T',-(H - W1*U1')*N1,'N',-N1);
    transfPairs(4)=struct('R',W2*U2','T',-(H - W2*U2')*N2,'N',-N2);
    nFrontPoints=zeros(1,4);
    for i=1:4
        transfPairs(i).lambda=epipolarTriangulateDepths(transfPairs(i).R,transfPairs(i).T,x1,x2);
        nFrontPoints(i) = sum(transfPairs(i).lambda(1,:)>0 & transfPairs(i).lambda(2,:)>0);
    end
    [~,idxMax]=max(nFrontPoints);

    N=transfPairs(idxMax).N;
    R=transfPairs(idxMax).R;
    T=transfPairs(idxMax).T;
    lambda=transfPairs(idxMax).lambda;
end
