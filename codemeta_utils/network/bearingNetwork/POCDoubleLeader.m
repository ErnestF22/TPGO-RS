function POCDoubleLeader
%resetRands();

A=zeros(3);
A(1,3)=1;
A(2,3)=1;
A=A+A';
t_node=testNetworkCreateStruct(A);

xTruth=[0 -1; 0 1; 1 1]';
t_node=bearingNetworkAddGroundTruth(t_node,xTruth);

xInit=xTruth;
xInit(:,3)=2*randn(2,1);
t_node.Ti=xInit;

%optsControl={};
optsControl={'leader',[1 2],zeros(2)};
funsBearings=bearingCostFunctions('cosine');

figure(1)
[t,x]=bearingNetworkEvolve(t_node,'tFinal',10,'funsBearings',funsBearings,...
    'optsControl',optsControl);
xFinal=x(:,:,end);
t_node.Ti=xFinal;

figure(2)
bearingNetworkPlot(t_node)
hold on
plotPoints(xInit,'rx')
plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
hold off
