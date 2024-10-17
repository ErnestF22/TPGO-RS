function POCHomographyFlowSegmentation
resetRands()
[X,G,NVec,idxX]=homFlowDatasetStructure(2);
NX=size(X,2);

flagDisplayScene=true;
if flagDisplayScene
    figure(1)
    plotGroups(X,idxX)
    hold on
    draw3dcameraFromG(G,'references')
    hold off
    axis equal
end

[x,dx,v,w]=homFlowDatasetFlow(X,G);
sigmaNoise=0.00;
dx=dx+sigmaNoise*randn(size(dx));

p=homogeneous(x);

flagDisplayField=true;
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
disp(norm(dx(:)-dxB(:),Inf))

%Build inverse of depths
lInv=zeros(size(l));
for iX=1:NX
    n=-nVecinG(:,idxX(iX));
    lInv(iX)=n'*p(:,iX);
end
disp(norm(lInv-l.^-1,Inf))

%Second way to build dx
dxC=zeros(size(dx));
for iX=1:NX
    xi=x(1,iX);
    yi=x(2,iX);
    alpha=homFlowParametersFromMotion(v,w,NVecinG(:,idxX(iX)));
    %    1  2  3  4  5  6     7     8
    Ai=[ 1 xi yi  0  0  0 xi*yi  xi^2;
         0  0  0  1 xi yi  yi^2 xi*yi;
         ];
        
    dxC(:,iX)=Ai*alpha;
end
disp(norm(dx(:)-dxC(:),Inf))

idxList=nearestNeighborGraphFromPoints(x,7);

alphaLocal=zeros(8,NX);
for iX=1:NX
    alphaLocal(:,iX)=homFlowParametersEstimate7pt(x(:,idxList(iX,:)),dx(:,idxList(iX,:)));
end
DAlpha=euclideanDistMatrix(alphaLocal,alphaLocal);
phi=@(dSq) exp(-0.01*dSq/2);
[idxAlpha,info]=quickshift_cluster(DAlpha,'threshold',0.01,'phi',phi);

figure(3)
nearestNeighborGraphPlot(x,idxList,[],'color',0.8*ones(1,3))
hold on
plotGroups(x,idxAlpha)
hold off

function overparametrization
dxC=zeros(size(dx));
thetaTruth=[-vec(v*NVecinG(1:3,1)');w];
for iX=1:NX
    M=[1 0 -x(1,iX); 0 1 -x(2,iX)];
    %lInv(iX)=n'*p(:,iX);
    %dxC(:,iX)=-M*v*n'*p(:,iX)+M*hat(p(:,iX))*w;
    Ai=[-kron(p(:,iX)',M) M*hat(p(:,iX))];
    n=-NVecinG(1:3,idxX(iX));
    dxC(:,iX)=Ai*[vec(v*n');w];
    %disp(norm(dx(:,iX)-dxC(:,iX),Inf))
end
%disp(norm(dx(:)-dxC(:),Inf))


A=zeros(2*NX,12);
idxA=reshape(1:2*NX,2,NX);
for iX=1:NX
    M=[1 0 -x(1,iX); 0 1 -x(2,iX)];
    Ai=[-kron(p(:,iX)',M) M*hat(p(:,iX))];
    A(idxA(:,iX),:)=Ai;
end
disp(norm(dx(:)-A*thetaTruth,Inf))


theta=A\dx(:);

% wEst=theta(10:12);
% vnEst=reshape(theta(1:9),3,3);
% [U,S,V]=svd(vnEst);
% vnEst=U(:,1)*S(1,1)*V(:,1)';

%disp([w wEst])
%disp([v*NVecinG(1:3,1)' vnEst])

keyboard

function a=vec(A)
a=A(:);
