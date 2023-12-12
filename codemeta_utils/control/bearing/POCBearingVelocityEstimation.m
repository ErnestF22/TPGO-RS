function POCBearingVelocityEstimation
resetRands()
costNameRanges='squared';
funsRanges=bearingCostFunctions(costNameRanges);
NNodes=4;

t_node=bearingNetworkBuildTestNetwork(NNodes);
EBearings=t_node.E;
ERanges=[1 2; 2 3];
NEdgesBearings=size(EBearings,1);
NEdgesRanges=size(ERanges,1);

xg=t_node.Titruth;
[ygBearings,nygBearings]=bearingNetworkComputeBearings(xg,EBearings);
[ygRanges,nygRanges]=bearingNetworkComputeBearings(xg,ERanges);

[x,~,~,dx]=real_randGeodFun(t_node.Ti);

nyBearings=@(t) bearingNetworkComputeRanges(x(t),EBearings);
yBearings=@(t) bearingNetworkComputeBearings(x(t),EBearings);
dyBearings=@(t) bearingNetworkComputeBearingsDerivative(dx,yBearings(t),nyBearings(t),EBearings);
nyRanges=@(t) bearingNetworkComputeRanges(x(t),ERanges);
yRanges=@(t) bearingNetworkComputeBearings(x(t),ERanges);

% A=BBearings'*inv(R)*BBearings;
v1=kron(ones(NNodes,1),[1;0])/sqrt(4);
v2=kron(ones(NNodes,1),[0;1])/sqrt(4);
P1=eye(2*NNodes)-v1*v1';
P2=eye(2*NNodes)-v2*v2';
P12=P1*P2*P1;
% Pdx=A*pinv(A);

% check_der(yBearings,@(t) bearingNetworkComputeBearingsDerivative(dx,yBearings(t),nyBearings(t),EBearings));

q=@(t) bearingNetworkComputeRangeResiduals(yRanges(t),ygRanges,nyRanges(t));
dq=@(t) bearingNetworkComputeRangeResidualsDerivatives(dx,yRanges(t),ygRanges,ERanges);

% check_der(q,dq)

d=size(dx,1);
R=@(t) bearingBuildR(nyBearings(t),d);
BBearings=@(t) bearingBuildBBearings(yBearings(t),NNodes,EBearings);
BRanges=@(t) bearingBuildBRanges(yRanges(t),ygRanges,NNodes,ERanges);
B=@(t) [BBearings(t);BRanges(t)];
D=@(t) blkdiag(inv(R(t)),eye(NEdgesRanges));
z=@(t) [reshape(yBearings(t),[],1); q(t)];
dz=@(t) D(t)*B(t)*dx(:);
%check_der(z,dz)

A=@(t) B(t)'*D(t)*B(t);

%dxRec=@(t) pinv(A(t))*B(t)'*dz(t);
%dxRec=@(t) pinv(A(t))*(BBearings(t)'*reshape(dyBearings(t),[],1)+BRanges(t)'*dq(t));

%uBearings=@(t) BBearings(t)'*reshape(dyBearings(t),[],1);
uBearings=@(t) reshape(controlBearings(dyBearings(t),NNodes,EBearings),[],1);
%uRanges=@(t) BRanges(t)'*dq(t);
uRanges=@(t) reshape(controlRanges(dq(t),ygRanges,NNodes,ERanges),[],1);
dxRec=@(t) pinv(A(t))*P12*(uBearings(t)+uRanges(t));

figure(1)
t=linspace(0,10,100);
plotfun(dxRec,t,'b')
hold on
plotfun(@(t) P12*dx(:),t,'rx')
hold off

disp([pinv(A(0))*A(0)-P12])

function uBearings=controlBearings(dy,NNodes,E)
d=size(dy,1);
uBearings=zeros(d,NNodes);

for iNode=1:NNodes
    uBearings(:,iNode)=-sum(dy(:,E(:,1)==iNode),2);
    uBearings(:,iNode)=uBearings(:,iNode)+sum(dy(:,E(:,2)==iNode),2);
end

function uRanges=controlRanges(dq,yg,NNodes,E)
d=size(yg,1);
uRanges=zeros(d,NNodes);
for iNode=1:NNodes
    iEdgeIn=E(:,1)==iNode;
    iEdgeOut=E(:,2)==iNode;
    uRanges(:,iNode)=-sum(yg(:,iEdgeIn).*(ones(d,1)*dq(iEdgeIn)'),2);
    uRanges(:,iNode)=uRanges(:,iNode)+sum(yg(:,iEdgeOut).*(ones(d,1)*dq(iEdgeOut)'),2);
end

function q=bearingNetworkComputeRangeResiduals(y,yg,ny)
q=ny.*bearingNetworkComputeBearingsCosines(y,yg)';

function dq=bearingNetworkComputeRangeResidualsDerivatives(dx,y,yg,E)
NNodes=size(dx,2);
Bq=bearingBuildBRanges(y,yg,NNodes,E);
dq=Bq*dx(:);

function Bq=bearingBuildBRanges(y,yg,NNodes,E)
d=size(y,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
Bq=zeros(NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);

    Bq(iEdge,idxNodes(:,iNode))=-yg(:,iEdge);
    Bq(iEdge,idxNodes(:,jNode))=yg(:,iEdge);
end

function dy=bearingDerivative(v,y,ny,E)
[d,NNodes]=size(v);
R=bearingBuildR(ny,d);
B=bearingBuildBBearings(y,NNodes,E);
dy=reshape(R\B*v(:),2,[]);

function R=bearingBuildR(nYij,d)
R=kron(diag(nYij),eye(d));

function B=bearingBuildBBearings(Yij,NNodes,E)
d=size(Yij,1);
NEdges=size(E,1);
idxNodes=reshape(1:d*NNodes,d,NNodes);
idxEdges=reshape(1:d*NEdges,d,NEdges);
B=zeros(d*NEdges,d*NNodes);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    
    PYij=eye(2)-Yij(:,iEdge)*Yij(:,iEdge)';
    B(idxEdges(:,iEdge),idxNodes(:,iNode))=-PYij;
    B(idxEdges(:,iEdge),idxNodes(:,jNode))=PYij;
end
