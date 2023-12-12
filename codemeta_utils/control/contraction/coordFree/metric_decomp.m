function [U,betap]=metric_decomp(M,alpha,beta)
% Use Schur Complement to convert [beta alpha'; alpha M] into its LDU
% decomposition
d=size(M,1);
betap=beta-alpha'*(M\alpha);
U=[1 zeros(1,d);-M\alpha eye(d)];