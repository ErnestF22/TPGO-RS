function [D] = eulBundle_metric(Z,X,Y,m,n)
% Compute the metric on TR^n
% INPUTS:
%   Z := Point of evaluation (p,u) given as a numeric matrix [6 x 1] (not
%       used)
%   X(p,u), Y(p,u) := vector fields on TR^n evaluated at (p,u)
%   m := an array of 3 values for the non natural metric on TR^n
%   n := scalar representing the dimension of the manifold
% OUTPUTS:
%   D := the nonNatural metric <X,Y> on TR^n

%% Extract useful parameters
m1 = m(1); m2 = m(2); m3 = m(3);
p = [eye(n) zeros(n)]*Z;
u = [zeros(n) eye(n)]*Z;
Xh = [eye(n) zeros(n)]*X;
Xv = [zeros(n) eye(n)]*X;
Yh = [eye(n) zeros(n)]*Y;
Yv = [zeros(n) eye(n)]*Y;

%% Compute the nonNatural metric
metric_Xh_Yh = m1*eul_metric(Xh,Yh);
metric_Xh_Yv = m2*eul_metric(Xh,Yv);
metric_Xv_Yh = m2*eul_metric(Xv,Yh);
metric_Xv_Yv = m3*eul_metric(Xv,Yv);

D = metric_Xh_Yh + metric_Xh_Yv + metric_Xv_Yh + metric_Xv_Yv;
end

