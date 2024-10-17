function POCFactorizationBearingLocalization_analysis
resetRands()
global AMeasurements
global bMeasurements
load('POCFactorizationBearingLocalization_data.mat')
disp('Original data')
residuals(W)
disp('Estimated matrix')
residuals(WEstimated)
disp('Uniform scaling of translation columns')
QTScaling=diag([10*rand*ones(1,nbNodes) ones(1,3)]);
residuals(W*QTScaling)

s=estimateScale(W,WEstimated);
QTScaling=diag([s*ones(1,nbNodes) ones(1,3)]);
WEstimated=WEstimated*QTScaling;
scatter(vec(reduceW(W)),vec(reduceW(WEstimated)),'.')
hold on
plot([-20 20],[-20,20])
hold off
disp([reduceW(W) reduceW(WEstimated)])
[U,S,V]=svdRank3(W);
[UEstimated,SEstimated,VEstimated]=svdRank3(WEstimated);
disp('U- and V-subpace angles between W and WEstimated')
disp(subspace(U,UEstimated))
disp(subspace(V,VEstimated))

%WTest=[U(:,1:2) U(:,3)]*diag([S(1,1),S(2,2),S(3,3)])*[V(:,1:2) V(:,3)]';
%residuals(WTest)
keyboard

function n=residuals(W)
global AMeasurements
global bMeasurements
disp(norm(AMeasurements*vec(W)-bMeasurements))

function s=estimateScale(W,WEstimated)
nbNodes=size(W,1)/3;
[~,~,V]=svd([vec(W(:,1:nbNodes)),vec(WEstimated(:,1:nbNodes))],'econ');
p=V(:,end);
s=-p(2)/p(1);

function WReduced=reduceW(W)
nbNodes=size(W,1)/3;
WReduced=W(:,1:nbNodes);

function [U,S,V]=svdRank3(W)
[U,S,V]=svd(W,'econ');
U=U(:,1:3);
S=S(1:3,1:3);
V=V(:,1:3);

