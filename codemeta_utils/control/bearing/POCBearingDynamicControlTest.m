function POCBearingDynamicControlTest

resetRands()
lambda=0;

alpha=[1 5];
alphak=1;

costNameBearings='cosine';
costNameRanges='squared';

funsBearings=bearingCostFunctions(costNameBearings);
funsRanges=bearingCostFunctions(costNameRanges);

t_node=bearingNetworkBuildTestNetwork();
NNodes=t_node.NNodes;

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

ERanges=t_node.Er;
ygRanges=t_node.Yrijtruth;
nygRanges=t_node.nYrijtruth;

xg=t_node.Titruth;
offset=[0;0];
x0=6*randn(2,t_node.NNodes);
x0=x0-(mean(x0,2)+offset)*ones(1,t_node.NNodes);

z0=[zeros(size(x0));x0];

dzControl=@(t,z)  closedLoop(z,EBearings,ygBearings,funsBearings,...
    ERanges,ygRanges,nygRanges,funsRanges,...
    lambda,alpha,alphak);
cEval=@(z) cost(z,EBearings,ygBearings,funsBearings,...
    ERanges,ygRanges,nygRanges,funsRanges,...
    alpha);


u=randn(2,NNodes);
z=@(t) [2*t*u;t^2*u];
t=linspace(-1,1,50);

switch 2
    case 1
        %check_der(@(t) costAndDer(cEval,z(t),dz(t,z(t))))
        figure(1)
        subplot(2,1,1)
        dz=@(t) [2*u; 2*t*u];
        check_der(@(t) costAndDer(cEval,z(t),dz(t)),'function',t)
        title('Function and derivative')
        subplot(2,1,2)
        plotfun(@(t) der(cEval,z(t),dzControl(t,z(t))),t)
        title('Inner product of control with gradient')
    case 2
        t1=t(10);
        z1=z(t1);
        dz1=reshape(dzControl(t1,z1),4,[]);
        [c1,gradc1]=cEval(z1);

        disp(gradc1(:)'*dz1(:))
        
        zGeod=@(t) z1+t*dz1+0.5*t^2*[zeros(2,NNodes);dz1(1:2,:)];
        figure(1)
        plotfun(@(t) cEval(zGeod(t)),linspace(0,1e-1,50));
end

%check_der(@(t) costAndDerPot(z(t),dz(t),EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)); 
% 
% function [c,dc]=costAndDerPot(z,dz,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)
% x=z(3:4,:);
% x2=z(1:2,:);
% dx=dz(3:4,:);
% dx2=dz(1:2,:);
% [yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
% [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
% 
% [c,gradPot]=bearingNetworkCostCombined(...
%     EBearings,ERanges,...
%     yBearings,yRanges,ygBearings,ygRanges,...
%     nyBearings,nyRanges,nygRanges,...
%     funsBearings,funsRanges,alpha);
% 
% c=c+0.5*x2(:)'*x2(:);
% 
% dc=[x2(:);gradPot(:)]'*[dx2(:);dx(:)];

function dc=der(cEval,z,dz)
[~,dc]=costAndDer(cEval,z,dz);

function [c,dc]=costAndDer(cEval,z,dz)
dz=reshape(dz,4,[]);
[c,gradc]=cEval(z);
dc=gradc(:)'*dz(:);

function [c,gradc]=cost(z,EBearings,ygBearings,funsBearings,ERanges,ygRanges,nygRanges,funsRanges,alpha)
x=z(3:4,:);

[yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
[yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);

[cPot,gradPot]=bearingNetworkCostCombined(...
    EBearings,ERanges,...
    yBearings,yRanges,ygBearings,ygRanges,...
    nyBearings,nyRanges,nygRanges,...
    funsBearings,funsRanges,alpha);

% u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
%     ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);

dx=z(1:2,:);
% dyBearings=bearingNetworkDerivative(dx,yBearings,nyBearings,EBearings);
% dqRanges=bearingNetworkComputeRangeResidualsDerivatives(dx,yRanges,ygRanges,ERanges);

cKin=0.5*sum(dx(:).^2);
c=cPot+cKin;

% u2=bearingNetworkControlDynamic(EBearings,dyBearings,...
%     ERanges,dqRanges,ygRanges,alpha);
% 
% u1b=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
%     ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);
% u2b=bearingNetworkControlDynamic(EBearings,dyBearings,...
%     ERanges,dqRanges,ygRanges,alpha);
% 
% u=alphak(1)*u1b+alphak(2)*u2b;
% dc=dx(:)'*gradPot(:)+dx(:)'*u(:);
gradc=[dx;gradPot];
%gradc=[zeros(size(dx(:)));alphak(1)*gradPot(:)];

%state z:
%   first two rows:     dx
%   second two rows:    x
function dz=closedLoop(z,EBearings,ygBearings,funsBearings,...
    ERanges,ygRanges,nygRanges,funsRanges,lambda,alpha,alphak)
z=reshape(z,4,[]);
[yBearings,dyBearings]=measurementBearings(z,EBearings);
[yRanges,nyRanges,dqRanges]=measurementRanges(z,ERanges,ygRanges);
u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,...
    ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha,alphak);

dz=model(z,u,lambda); %dz=[u;z(1:2,:)];
%dz=[zeros(2,size(z,2)); u];
dz=dz(:);

function [y,dy]=measurementBearings(z,E)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dy=bearingNetworkDerivative(dx,y,ny,E);

function [y,ny,dq]=measurementRanges(z,E,yg)
x=z(3:4,:);
dx=z(1:2,:);
[y,ny]=bearingNetworkComputeBearings(x,E);
dq=bearingNetworkComputeRangeResidualsDerivatives(dx,y,yg,E);

function u=control(EBearings,yBearings,dyBearings,ygBearings,funsBearings,...
    ERanges,yRanges,dqRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha,alphak)
u1=bearingNetworkControlDirect(EBearings,yBearings,ygBearings,funsBearings,...
    ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges,alpha);
u2=bearingNetworkControlDynamic(EBearings,dyBearings,...
    ERanges,dqRanges,ygRanges,alpha);

u=u1+alphak*u2;

function dz=model(z,u,lambda)
m=1;
dz=zeros(size(z));
dz(1:2,:)=1/m*u-lambda*z(1:2,:);
dz(3:4,:)=z(1:2,:);

function dYij=bearingNetworkDerivative(v,Yij,nYij,E)
[d,NNodes]=size(v);
R=bearingBuildR(nYij,d);
B=bearingBuildB(Yij,NNodes,E);
dYij=reshape(R\B*v(:),2,[]);

function R=bearingBuildR(nYij,d)
R=kron(diag(nYij),eye(d));

function B=bearingBuildB(Yij,NNodes,E)
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
