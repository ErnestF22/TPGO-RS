%Solve the low-rank least squares problem with ADMM and SVD projection
%Solves the following problem: find a matrix X by 
%min ||A*vec(X)-b||^2 
%subject to rank(X)=k
% WARNING: THIS POC DOES NOT GIVE THE RIGHT SOLUTION!
function POCFixedRankLSwithADMMSVD

sz=[30,20];
k=3;
nbConstraints=7;

X=randn(sz(1),k)*randn(k,sz(2));

A=randn(nbConstraints,prod(sz));
b=A*vec(X);

[X,U,output]=ADMM(A,b,sz,k);

figure(1)
subplot(3,1,1)
semilogy(output.cost.linearConstraints);
subplot(3,1,2)
semilogy(output.cost.residuals);
subplot(3,1,3)
semilogy(output.svals');

svd(U)

function [X,U,output]=ADMM(A,b,C,D,sz,k)
%parameters
rho=1;
maxIt=100;

%elements of the quadratic part from the LS problem
P=A'*A;
q=A'*b;
I=eye(size(P));

%initialization
%X is the final solution
X=reshape(A\b,sz(1),sz(2));
%U contains the Lagrange multipliers
U=zeros(sz);
%Z contains the projection of X onto the fixed rank set
Z=projRank(X,k);

%initialization of debug info
output.cost.linearConstraints=zeros(1,maxIt);
output.cost.residuals=zeros(1,maxIt);
output.svals=zeros(min(sz),maxIt);

%ADMM iterations
for it=1:maxIt
    v=vec(Z-U);
    X=reshape((P+rho*I)\(rho*v-q),sz(1),sz(2));
    Z=projRank(X+U,3);
    U=U+X-Z;
    %update debug info
    output.cost.linearConstraints(it)=norm(A*vec(X)-b,2)^2;
    output.cost.residuals(it)=norm(X-Z,'fro');
    output.svals(:,it)=svd(X);
end

%project a matrix onto the set of rank(X)=k matrices
function Z=projRank(X,k)
[U,S,V]=svd(X);
Z=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';
