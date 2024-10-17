function POCHomographyFlowEstimation
resetRands()
[X,G,NVec,idxX]=homFlowDatasetStructure(1);
NX=size(X,2);

flagDisplayScene=false;
if flagDisplayScene
    figure(1)
    plotGroups(X,idxX)
    hold on
    draw3dcameraFromG(G,'references')
    hold off
    axis equal
end

[x,dx,v,w]=homFlowDatasetFlow(X,G);


p=homogeneous(x);

flagDisplayField=false;
if flagDisplayField
    figure(2)
    plotField(x,dx)
    axis equal
end

%First way to build dx
s=[ones(1,NX); zeros(1,NX);-x(1,:)];
r=[zeros(1,NX); ones(1,NX);-x(2,:)];
l=projectGetDepthsFromG(G,X,'references');
NVecinG=rigidTransformG(G,NVec,'references','planes','wc');
nVecinG=planeNVecToNScaled(NVecinG);

dxB=-[
    (v'*s)./l+sum(s.*(hat(w)*p));
    (v'*r)./l+sum(r.*(hat(w)*p));
    ];
%disp(norm(dx(:)-dxB(:),Inf))

%Build inverse of depths
lInv=zeros(size(l));
for iX=1:NX
    n=-nVecinG(:,idxX(iX));
    lInv(iX)=n'*p(:,iX);
end
%disp(norm(lInv-l.^-1,Inf))

%build matrix going from n*v' to alpha
VNToAlphaMat=zeros(8,9);
idxMat=reshape(1:9,3,3);

VNToAlphaMat(1,idxMat(3,1))=1;
VNToAlphaMat(2,idxMat(1,1))=1;
VNToAlphaMat(2,idxMat(3,3))=-1;
VNToAlphaMat(3,idxMat(2,1))=1;
VNToAlphaMat(4,idxMat(3,2))=1;
VNToAlphaMat(5,idxMat(1,2))=1;
VNToAlphaMat(6,idxMat(2,2))=1;
VNToAlphaMat(6,idxMat(3,3))=-1;
VNToAlphaMat(7,idxMat(2,3))=-1;
VNToAlphaMat(8,idxMat(1,3))=-1;


%Second way to build dx
dxC=zeros(size(dx));
n=nVecinG(:,idxX(iX));
nv=n*v';

alpha=homFlowParametersFromVN(v,NVecinG(:,idxX(iX)));
%disp([alpha VNToAlphaMat*vn(:)])

A=[];
b=[];
A1Fact=[];
A2Fact=[];
NXEst=4;
dxEst=dx(:,1:NXEst);
M=[];
for iX=1:NXEst
    xi=x(1,iX);
    yi=x(2,iX);
    %    1  2  3  4  5  6     7     8
%     Ai=[ 1 xi yi  0  0  0 xi*yi  xi^2;
%          0  0  0  1 xi yi  yi^2 xi*yi;
%          ]*VNToAlphaMat;
    A2i=[eye(2) -x(:,iX)]';
    A1i=[x(:,iX)' 1];
    Ai=kron(A2i',A1i);
    bi=[-w(2)+yi*w(3)+xi*yi*w(1)-xi^2*w(2);
        w(1)-xi*w(3)+yi^2*w(1)-xi*yi*w(2)];
    Mi=(dxEst(:,iX)-bi)';
    
    %disp(Ai*vec(vn)-vec(A1i*vn*A2i));
    %disp(A1i*vn*A2i-Mi)
    disp(A1i*nv*A2i-Mi)
    A=[A;Ai];
    b=[b;bi];
    M=[M;Mi];
    dxC(:,iX)=Ai*nv(:)+bi;
end
%disp(norm(dx(:)-dxC(:),Inf))
disp(norm(dxEst(:)-b-A*nv(:),Inf))

nvEst=reshape(pinv(A)*(dxEst(:)-b),3,3);
s=svd(nvEst);
nvEst=nvEst-s(2)*eye(3);
disp(norm(dxEst(:)-b-A*nvEst(:),Inf))
disp([nv(:) nvEst(:)])

%lambdaOpt=mean(diag(vn-vnEst));
% funPlot(@(l) costBestRank1Approx(vnEst,l),linspace(-1,1,200))
% hold on
% ax=axis();
% plot([lambdaOpt lambdaOpt], ax(3:4),'k--')
% hold off

function c=costBestRank1Approx(A,lambda)
c=[0 1 1]*svd(A+lambda*eye(3));



