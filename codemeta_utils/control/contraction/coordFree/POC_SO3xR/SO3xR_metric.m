function [D] = SO3xR_metric(Z,X,Y)
% Compute the natural metric on SO3xR. 
%%%%%Assume all matrices are numerics!%%%%%
% INPUTS:
%   Z := Point of evaluation (R,p)
%   X,Y := Vector fields
% OUTPUTS:
%   D := result of <X,Y> on SO3xR

% Extract components
Z_SO3 = extractComp(Z,1,3,1,3);
X_SO3 = extractComp(X,1,3,1,3);
Y_SO3 = extractComp(Y,1,3,1,3);
X_R = extractComp(X,4,4,1,1);
Y_R = extractComp(Y,4,4,1,1);

% Compute metric on SO3
D_SO3 = rot_metric(Z_SO3,X_SO3,Y_SO3);
% Compute metric on R
D_R = 1/2*X_R*Y_R;

% Sum metric
D = D_SO3+D_R;
end

