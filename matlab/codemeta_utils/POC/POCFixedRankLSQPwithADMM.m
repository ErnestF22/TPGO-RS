%Solve the low-rank least squares problem with ADMM and SVD projection
%Solves the following problem: find a matrix X by 
%min ||A*vec(X)-b||^2 
%subject to 
%   C*vec(X)>D
%   rank(X)=k
% WARNING: THIS POC DOES NOT GIVE THE RIGHT SOLUTION!
% Instead, use fixedRankLSQP.m
function POCFixedRankLSQPwithADMM
resetRands(0)
sz=[20,10];
k=3;
nbConstraints=7;
nbInequalities=3;

X=randn(sz(1),k)*randn(k,sz(2));

A=randn(nbConstraints,prod(sz));
b=A*vec(X);

C=randn(nbInequalities,prod(sz));
d=randn(nbInequalities,1);

[X,U,output]=ADMM(A,b,C,d,sz,k);

figure(1)
subplot(5,1,1)
semilogy(output.cost.linearConstraints);
title('LS cost')
subplot(5,1,2)
semilogy(output.cost.fixedRankResiduals);
title('Residual rank constraint')
subplot(5,1,3)
semilogy(output.cost.inequalitiesResiduals);
title('Residuals inequality constraints')
subplot(5,1,4)
semilogy(output.svals');
title('Svals of X')
subplot(5,1,5)
plot(output.cost.inequalitiesNbActive)
title('Nb. of active inequalities')
svd(U)

function [X,UFixedRank,output]=ADMM(A,b,C,d,sz,k)
%parameters
rho=1;
maxIt=100;
tolInequalities=1e-6;
nbInequalities=size(C,1);

%elements of the quadratic part from the LS problem
P=A'*A;
q=A'*b;
I=eye(size(P));

%normalize inequalities
[C,d]=linearInequalitiesNormalize(C,d);

%initialization
%X is the final solution
X=reshape(A\b,sz(1),sz(2));
%UFixedRank contains the Lagrange multipliers for the fixed rank constraint
UFixedRank=zeros(sz);
%UInequality contains the Lagrange multipliers for each inequality
%constraints
UInequalities=zeros(prod(sz),nbInequalities);
%Z contains the projection of X onto the fixed rank set
Z=fixedRankProject(X,k);
%Y contains the projection of X on each one of the inequalities
Y=linearInequalitiesProject(vec(X),C,d,'normalized');

%initialization of debug info
output.cost.linearConstraints=zeros(1,maxIt);
output.cost.fixedRankResiduals=zeros(1,maxIt);
output.cost.inequalitiesResiduals=zeros(nbInequalities,maxIt);
output.cost.inequalitiesNbActive=zeros(1,maxIt);
output.svals=zeros(min(sz),maxIt);

%ADMM iterations
for it=1:maxIt
    vFixedRank=vec(Z-UFixedRank);
    vInequalities=Y-UInequalities;
    X=reshape((P+(nbInequalities+1)*rho*I)\(rho*(vFixedRank+sum(vInequalities,2))-q),sz(1),sz(2));
    Z=fixedRankProject(X+UFixedRank,3);
    Y=linearInequalitiesProject(vec(X),C,d,'normalized');
    UFixedRank=UFixedRank+X-Z;
    UInequalities=UInequalities+vec(X)-Y;
    
    %update debug info
    output.cost.linearConstraints(it)=norm(A*vec(X)-b,2)^2;
    output.cost.fixedRankResiduals(it)=norm(X-Z,'fro');
    output.cost.inequalitiesResiduals(:,it)=min(C*vec(X)-d,0);
    output.cost.inequalitiesNbActive(:,it)=sum(C*vec(X)-d<tolInequalities);
    output.svals(:,it)=svd(X);
end

%project a matrix onto the set of rank(X)=k matrices
function Z=fixedRankProject(X,k)
[U,S,V]=svd(X);
Z=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';

