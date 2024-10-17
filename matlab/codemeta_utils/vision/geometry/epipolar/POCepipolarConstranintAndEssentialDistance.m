function POCepipolarConstranintAndEssentialDistance
%resetRands();
[G,X]=testNetworkCreateAbsolutePoses(3);
%X=X(:,1:5);
G=G(:,:,1:2);
%testNetworkDisplay(G,'points',X);
x=projectFromG(G,X);
%plot(x(1,:,1),x(2,:,1),'*')
QTruth=essential_fromG(G(:,:,1),G(:,:,2),'poses');
Q0=essential_randn(QTruth,0.01);

QEst=essential_minimizeEpipolarCost(Q0,x(:,:,1),x(:,:,2),...
    'optsMinimize',{'stopThreshold',1e-15,'displayIt'});

[eEst,~,HeEst]=essential_evaluateEpipolarCost(QEst,x(:,:,1),x(:,:,2));
[eTruth,~,HeTruth]=essential_evaluateEpipolarCost(QTruth,x(:,:,1),x(:,:,2));
disp('Distance between QTruth/Q0 and QTruth/QEst')
disp([essential_dist(QTruth,Q0) essential_dist(QTruth,QEst)])
disp('Epipolar cost evaluation QEst/QTruth')
disp([eEst;eTruth]);
subplot(1,2,1)
essential_disp(QTruth);
subplot(1,2,2)
essential_disp(QEst);
disp('Singular values of the Hessian of the total cost QEst/QTruth')
[svd(sum(HeEst,3)) svd(sum(HeTruth,3))]
