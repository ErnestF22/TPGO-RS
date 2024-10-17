function POCessential_optimalSfM
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
X8pt=Qcol2row(Q8pt);

problem.M=essentialfactory();
problem.cost=@(X) f(X,x);
problem.grad=@(X) gradf(X,x);
problem.hess=@(X,S) hessf(X,S,x);

% % Numerically check the differentials.
% checkgradient(problem); pause;
% checkhessian(problem); pause;

XOpt=trustregions(problem,X8pt);

QOpt=Qrow2col(XOpt);

disp('Sampson cost (8pt/Opt)')
disp([essential_evaluateEpipolarCostSampsonSq(Q8pt,x) essential_evaluateEpipolarCostSampsonSq(QOpt,x)])

disp('Essential distance from QTruth (8pt/Opt)')
disp([essential_dist(QTruth,Q8pt,'signed') essential_dist(QTruth,QOpt,'signed')])


function X=Qcol2row(Q)
X=[Q(1:3,:,:) Q(4:6,:,:)];

function Q=Qrow2col(X)
Q=[X(:,1:3,:); X(:,4:6,:)];

function c=f(X,x)
c=essential_evaluateEpipolarCostSampsonSq(Qrow2col(X),x);

function g=gradf(X,x)
[~,gQ]=essential_evaluateEpipolarCostSampsonSq(Qrow2col(X),x);
g=Qcol2row(gQ);

function h=hessf(X,S,x)
[~,~,hessOpQ]=essential_evaluateEpipolarCostSampsonSq(Qrow2col(X),x,'symmetricHess');
h=Qcol2row(hessOpQ(Qrow2col(S)));

