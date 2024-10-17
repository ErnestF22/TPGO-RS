%Decomposes R such that V*S*V'=R and S is of the form
%   [cos -sin 0]
%S= [sin cos  0]
%   [ 0    0  1]

%TODO: numerically this is not the best and it could be improved
function [V,S]=rot_realdecomp(R)
[V,S]=eig(R);
M=[-1/sqrt(2) -1/sqrt(2) 0;i/sqrt(2) -i/sqrt(2) 0; 0 0 1];
V=V*inv(M);
S=M*S*inv(M);
