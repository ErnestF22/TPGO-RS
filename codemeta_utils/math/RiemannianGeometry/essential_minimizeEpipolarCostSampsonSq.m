%Optimize Sampson's epipolar cost on the essential manifold
%function [Q,cost,info]=essential_minimizeEpipolarCostSampsonSq(x,QInit)
%Inputs
%   x       [2 x Nx x 2] array with the image points for the two views
%   QInit   [6 x 3 x NQ] Initial guess for the optimization. If multiple
%           points are specified, the algorithm is run multiple times and
%           multiple results are returned.
function [Q,cost,info]=essential_minimizeEpipolarCostSampsonSq(x,QInit)
if ~exist('QInit','var')
    QInit=[];
end

XInit=essential_col2row(QInit);

%setup the manopt problem
problem.M=essentialfactory();
problem.cost=@(X) f(X,x);
problem.grad=@(X) gradf(X,x);
problem.hess=@(X,S) hessf(X,S,x);

NQ=size(XInit,3);
X=zeros(size(XInit));
cost=zeros(1,NQ);
info=cell(1,NQ);

for iQ=1:NQ
    [X(:,:,iQ),cost(iQ),info{iQ}]=trustregions(problem,XInit(:,:,iQ),struct('verbosity',0));
end

Q=essential_row2col(X);

%wrapper functions for computing cost, gradient and Hessian as required by
%manopt
%These could be made much more efficient by storing values in the algorithm
function c=f(X,x)
c=essential_evaluateEpipolarCostSampsonSq(essential_row2col(X),x);

function g=gradf(X,x)
[~,gQ]=essential_evaluateEpipolarCostSampsonSq(essential_row2col(X),x);
g=essential_col2row(gQ);

function h=hessf(X,S,x)
[~,~,hessOpQ]=essential_evaluateEpipolarCostSampsonSq(essential_row2col(X),x,'symmetricHess');
h=essential_col2row(hessOpQ(essential_row2col(S)));

