function [ D ] = rotBundle_metric_sasaki( Z, X, Y )
% Compute the metric on the tangent bundle TSO(n). The induced Sasaki
% metric on TSO(n) is given on pg 14-15 in Gudmundsson and Kappos.
% INPUTS:
%   Z := Point of evaluation (p,u) given as a numeric matrix (IE [6 x 3]
%       for TSO(3)
%   X, Y := Vector fields on TSO(n) given as numeric matrices (IE: [6 x 3]
%       for TSO(3)
% OUTPUTS:
%   D := Result of <X,Y> on TSO(n) 

%% Split X and Y into vertical and horizontal components
[R, ~] = rotBundle_split(Z);
Xh = rotBundle_extractHoriz(Z, X);
Xv = rotBundle_extractVert(Z, X);
Yh = rotBundle_extractHoriz(Z, Y);
Yv = rotBundle_extractVert(Z, Y);

%% <Xh, Yh> on TSO(n) == <X,Y> on SO(n)
metric_Xh_Yh = rot_metric(R, Xh, Yh);

%% <Xv, Yh> on TSO(n) == 0
metric_Xv_Yh = 0;

%% <Xv, Yv> on TSO(n) == <X,Y> on SO(n)
metric_Xv_Yv = rot_metric(R, Xv, Yv);

%% Resulting D is the sum of the three components
D = metric_Xh_Yh + metric_Xv_Yh + metric_Xv_Yv;
end

