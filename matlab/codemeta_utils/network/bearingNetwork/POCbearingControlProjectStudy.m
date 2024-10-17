function POCbearingControlProjectStudy
%try to find local minima for the cost with projection
resetRands()


%make the network
A=adjgallery(4,'full');
t_node=testNetworkCreateStruct(A);

xtruth=[1 0; 0 0; -1 0.1; -1 -0.1]';
t_node=bearingNetworkAddGroundTruth(t_node,xtruth);


xInit=xtruth;
thetaInit=pi/4;
xInit(:,1)=[cos(thetaInit); sin(thetaInit)];
%xInit(:,[3 4])=[2.1 -0.1; 2 0.1]'; 
t_node.Ti=xInit;

optsControl={};
%optsControl={'leader',2,[0;0]};
%optsControl={'project',1};

% A=randn(2,2,t_node.NNodes);
% T=multiprod(multitransp(A),A);
% for iT=1:size(T,3)
%     T(:,:,iT)=T(:,:,iT)/det(T(:,:,iT));
% end
% optsControl={'preconditioner',T};

figure(1)
[t,x]=bearingNetworkEvolve(t_node,'tFinal',50,...
    'optsControl',optsControl);
xFinal=x(:,:,end);
t_node.Ti=xFinal;

figure(2)
bearingNetworkPlot(t_node)
hold on
plotPoints(xInit,'rx')
plot(squeeze(x(1,:,:))',squeeze(x(2,:,:))','-.','color',[1 0.75 0])
hold off
