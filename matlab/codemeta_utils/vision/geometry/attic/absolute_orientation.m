%function [R,t]=absolute_orientation(x1,x2)
%Solve the absolute orientation problem. Given two sets of corresponding 3D
%points x1 and x2, find the best R and t which minimize ||x2-R*x1-t||^2
function [R,t]=absolute_orientation(x1,x2)
n=size(x1,2);

x1bar=mean(x1,2);
x2bar=mean(x2,2);

x11=x1-x1bar*ones(1,n);
x21=x2-x2bar*ones(1,n);

M=x11*x21';

[U,S,V]=svd(M);
R=V*U';

t=x2bar-R*x1bar;
