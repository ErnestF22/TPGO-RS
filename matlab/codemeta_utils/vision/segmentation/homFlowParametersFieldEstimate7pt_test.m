function homFlowParametersFieldEstimate7pt_test
resetRands()
[X,G,~,idxX]=homFlowDatasetStructure(2,'NPointsPlane',80);
[x,dx]=homFlowDatasetFlow(X,G,'w',zeros(3,1));
sigmaNoise=0.03;
dx=dx+sigmaNoise*randn(size(dx));
[idxList,kList]=nearestNeighborGraphFromPoints(x,10,'symmetric');

lambdaBase=[5 1 1];
lambda2=10.^[-5 -3 0 2];
tic
[alphaHat,alpha,lambdaSchedule,alphaLocal]=homFlowParametersFieldEstimate7pt(x,dx,idxList,kList,...
    'lambda',lambdaBase,'lambda2',lambda2,'collect','maxItAlternation',10);
toc

E=homFlowParametersFieldEstimate7ptEnergy(alphaHat,alpha,x,dx,idxList,kList,...
    'lambdaSchedule',lambdaSchedule);

figure(1)
semilogy(E)
disp('Final energy')
disp(E(end))

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
