function essential_minimizeEpipolarCostSampsonSq_test
[G,X]=testNetworkCreateAbsolutePoses(3);
%X=X(:,1:5);
G=G(:,:,1:2);
testNetworkDisplay(G,'points',X);
x=projectFromG(G,X);

sigma=0.03;
x=x+sigma*randn(size(x));
x1=x(:,:,1);
x2=x(:,:,2);

QTruth=essential_fromG(G(:,:,1),G(:,:,2),'poses');

E8pt=epipolarEstimateE8pt(x1,x2);
Q8pt=essential_solveFlipAmbiguity(essential_fromE(E8pt),x);

tic
QOpt=essential_minimizeEpipolarCostSampsonSq(x,Q8pt);
toc

disp('Sampson cost (8pt/Opt)')
disp([essential_evaluateEpipolarCostSampsonSq(Q8pt,x) essential_evaluateEpipolarCostSampsonSq(QOpt,x)])

disp('Essential distance from QTruth (8pt/Opt)')
disp([essential_dist(QTruth,Q8pt,'signed') essential_dist(QTruth,QOpt,'signed')])
