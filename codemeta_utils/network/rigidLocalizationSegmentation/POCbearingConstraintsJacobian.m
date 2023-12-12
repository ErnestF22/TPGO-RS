function POCbearingConstraintsJacobian
resetRands(2);

graphNum=2;
dimData=2;
N=5;

% typeConstraints='relativeBearings';
% typeConstraints='fixedBearings';
typeConstraints='distances';

x=randn(dimData,N);
A=zeros(N);

switch graphNum
    case 1 % two separate components
        A(1,2)=1;
        A(3,4)=1;
        A(3,5)=1;
        A(4,5)=1;
    case 2 % two separate 2D rigid components
        A(1,2)=1;
        A(1,3)=1;
        A(2,3)=1;
        A(3,4)=1;
        A(3,5)=1;
        A(4,5)=1;
    case 3 % one 2D rigid component
        A(1,2)=1;
        A(1,3)=1;
        A(2,3)=1;
        A(3,4)=1;
        A(3,5)=1;
        A(4,5)=1;
        A(1,4)=1;
    case 4 % two 3D rigid component
        A(1,2)=1;
        A(1,3)=1;
        A(2,3)=1;
        A(3,4)=1;
        A(3,5)=1;
        A(4,5)=1;
        A(1,4)=1;
        A(1,5)=1;
    otherwise
        error('graphNum not valid')
end
A=A+A';

[C,E]=adj2incmatrix(A);%,'undirected');

switch typeConstraints
    case 'fixedBearings'
        dimConstraints=dimData;
        dimCoordinates=dimData;
    case 'distances'
        dimConstraints=1;
        dimCoordinates=dimData;        
    case 'relativeBearings'
        dimConstraints=dimData;
        switch dimData
            case 2
                dimCoordinates=3;
            case 3
                dimCoordinates=6;
        end
        R=rot_randn(eye(dimData),100,N);%repmat(eye(2),[1 1 N]);
        t=generateRelativeBearingMeasurements(x,R,E);
        %t=t+0.1*randn(size(t));
end


figure(1)
graphShow(E,x)
% axis([-1.5 1.5 -1.5 1.5])
% axis([0 6 0 6])


NEdges=size(E,1);
J=[];
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    xi=x(:,iNode);
    xj=x(:,jNode);
    switch typeConstraints
        case 'fixedBearings'
            Jr=blockRowJacobianFixedBearingsConstraints(N,iNode,jNode,xi,xj);
        case 'distances'
            Jr=blockRowJacobianDistanceConstraints(N,iNode,jNode,xi,xj);
        case 'relativeBearings'
            Ri=R(:,:,iNode);
            tij=t(:,iEdge);
            Jr=blockRowJacobianRelativeBearingsConstraints(N,iNode,jNode,xi,xj,Ri,tij);
        otherwise
            error('typeConstraints not valid');
    end
    J=[J; Jr];
end

disp('size(J)')
disp(size(J))
disp('rank(J)')
disp(rank(J))
disp('bN=null(J)')
bN=null(J);
displayMatrixInRowBlocks(bN,dimCoordinates)
centerNode=3;
disp(['bNCentered at node ' num2str(centerNode) ' and row-normalized'])
bNCentered=centerNullBasis(bN,centerNode,dimCoordinates);
bNCenterdNormalizd=cnormalize(bNCentered')';
displayMatrixInRowBlocks(bNCenterdNormalizd,dimCoordinates)

% disp('J*common translation')
% bT=cnormalize(kron(ones(N,1),eye(2)));
% disp(J*bT)
% disp('Non-trivial null(J)')
% [U,S,V]=svd((eye(N*2)-bT*bT')*bN);
% bN2=U(:,1:rank(S));
% % % bN2=cnormalize(bN2')';
% % disp(bN2)
% % figure(2)
% % plot(bN2(:,1),bN2(:,2),'*','MarkerSize',15)
% % disp('Indicator on all squared distances for embedded points')
% % allDists=computeAllEuclideanDistancesSq(bN2');
% % indAllDists=allDists<1e-6;
% % disp(indAllDists)
% % A1=indAllDists(1:2:end,1:2:end);
% % A2=indAllDists(2:2:end,2:2:end);
% % disp('Embedding adj matrices')
% % disp(A1)
% % disp(A2)
% disp('bN2Centered at node 3 and row-normalized')
% bN2Centered=centerNullBasis(bN2,3);
% bN2CenterdNormalizd=cnormalize(bN2Centered')';
% disp(bN2CenterdNormalizd)

% keyboard

save([mfilename '_data'])

function bNCentered=centerNullBasis(bN,iNode,dimCoordinates)
%subtract iNode-th block from all the others
if ~exist('dimCoordinates','var')
    dimCoordinates=2; %dimension of the coordinates for each node (i.e., of the blocks)
end
bNCentered=bN;
NNodes=size(bN,1)/dimCoordinates;
NiNode=bN((dimCoordinates*(iNode-1)+1):(dimCoordinates*iNode),:);
for jNode=1:NNodes
    bNCentered((dimCoordinates*(jNode-1)+1):(dimCoordinates*jNode),:)=bN((dimCoordinates*(jNode-1)+1):(dimCoordinates*jNode),:)-NiNode;
end

function Jr=blockRowJacobianFixedBearingsConstraints(N,iNode,jNode,xi,xj)
dimData=length(xi);
Jr=zeros(dimData,dimData*N);
J=singleJacobianBearingsConstraintsTranslation(xi,xj);
Jr(:,(dimData*(iNode-1)+1):dimData*iNode)=J;
Jr(:,(dimData*(jNode-1)+1):dimData*jNode)=-J;

function Jr=blockRowJacobianDistanceConstraints(N,iNode,jNode,xi,xj)
dimData=length(xi);
Jr=zeros(1,dimData*N);
J=(xi-xj)'/norm(xi-xj);
Jr(:,(dimData*(iNode-1)+1):dimData*iNode)=J;
Jr(:,(dimData*(jNode-1)+1):dimData*jNode)=-J;

function Jr=blockRowJacobianRelativeBearingsConstraints(N,iNode,jNode,xi,xj,Ri,tij)
dimData=length(xi);
switch dimData
    case 2
        dimCoords=3;
    case 3
        dimCoords=6;
end
Jr=zeros(dimData,dimCoords*N);
JTranslation=singleJacobianBearingsConstraintsTranslation(xi,xj);
switch dimData
    case 2
        JRotation=Ri*[0 -1; 1 0]*tij;
    case 3
        JRotation=-Ri*hat(tij);
end
J=[JRotation JTranslation];
Jr(:,(dimCoords*(iNode-1)+1):dimCoords*iNode)=J;
Jr(:,(dimCoords*(jNode-1)+1):dimCoords*jNode)=-J;

function J=singleJacobianBearingsConstraintsTranslation(xi,xj)
dimData=length(xi);
n=norm(xi-xj);
J=(eye(dimData)-(xi-xj)*(xi-xj)'/n^2)/n;

function t=generateRelativeBearingMeasurements(x,R,E)
dimData=size(x,1);
NEdges=size(E,1);
t=zeros(dimData,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    xi=x(:,iNode);
    xj=x(:,jNode);
    dij=xi-xj;
    Ri=R(:,:,iNode);
    t(:,iEdge)=Ri'*dij/norm(dij);
end

