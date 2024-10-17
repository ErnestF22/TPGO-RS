function homFlowParametersFieldEstimate4ptTranslation_test
resetRands()
[X,G,NVec,idxX]=homFlowDatasetStructure(1,'NPointsPlane',80);
w=zeros(3,1);
[x,dx,v]=homFlowDatasetFlow(X,G,'w',w);
NVecinG=rigidTransformG(G,NVec,'references','planes','wc');
n=planeNVecToNScaled(NVecinG);


sigmaNoise=0.03;
dx=dx+sigmaNoise*randn(size(dx));
[idxList,kList]=nearestNeighborGraphFromPoints(x,20,'symmetric');

lambdaBase=[5 1 1];
lambda2=10.^[-5 -3 0 2];
tic
[alphaHat,alpha,lambdaSchedule,alphaLocal]=homFlowParametersFieldEstimate4ptTranslation(x,dx,w,idxList,kList,...
    'lambda',lambdaBase,'lambda2',lambda2,'collect','maxItAlternation',10);
toc

% E=homFlowParametersFieldEstimate7ptEnergy(alphaHat,alpha,x,dx,idxList,kList,...
%     'lambdaSchedule',lambdaSchedule);
% 
% figure(1)
% semilogy(E)
% disp('Final energy')
% disp(E(end))

figure(2)
nearestNeighborGraphPlot(x,idxList,kList,'color',0.8*ones(1,3))
hold on
plotGroups(x,idxX)
hold off

flagDisplayField=true;
if flagDisplayField
    figure(3)
    plotField(x,dx)
    axis equal
end

figure(4)
subplot(2,1,1)
plot(alphaLocal')
subplot(2,1,2)
plot(alphaHat(:,:,end)')
hold on
plot(repmat([n;v],1,size(alphaHat,2)),'k')
hold off

save([mfilename '_data'])