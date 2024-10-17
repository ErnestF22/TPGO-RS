function POCbcnFranchi_formation
resetRands(2)
%%
NNodes=5;
d=2;
t_node.Titruth=randn(d,NNodes);
t_node.Ti=[t_node.Titruth(:,1:2) randn(d,NNodes-2)];
[t_node.Yijtruth,t_node.Gammaitruth,t_node.E]=bnFranchi_computeBearingAndGammaFormation(t_node.Titruth);
handle_ode=@(t,x) controller(t,x,t_node);
[t,x]=ode15s(handle_ode,[0 100],t_node.Ti(:));
x=reshape(x',d,NNodes,[]);
t_node.Ti=x(:,:,end);
bearingNetworkPlot(t_node)
hold on
bearingNetworkPlotTrajectories(x)
hold off


%%

function dx=controller(t,x,t_node)
yd=t_node.Yijtruth;
d=size(yd,1);
gd=t_node.Gammaitruth;
E=t_node.E;
k=1;
x=reshape(x,d,[]);
NNodes=size(x,2);
[y,g]=bnFranchi_computeBearingAndGammaFormation(x);
y12=y(:,1);
y12d=yd(:,1);
dx=zeros(size(x));
dx(:,2)=k*orthComplementProjector(y12)*y12d;
for iNode=3:NNodes
    idxy=2+2*iNode-5;
    y1i=y(:,idxy);
    y1id=yd(:,idxy);
    dx(:,iNode)=k*(gd(iNode-2)*y1id-g(iNode-2)*y1i);
end
dx=reshape(dx,[],1);

