function POCSymmBlockSkew
%Study of bounds on evals of a matrix of the form
% [Z A; A' Z], where A is skew symmetric

a=randn(3,1);
A=hat3(a);
Z=zeros(3);
M=[Z A; A' Z];
disp(max(eig(M-norm(a)*eye(6))))
