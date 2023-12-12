function [ D ] = rotBundle_metric_nonNatural( Z, X, Y, m )
% Compute the metric on the tangent bundle TSO(n). The induced Sasaki
% metric on TSO(n) is given on pg 14-15 in Gudmundsson and Kappos.
% INPUTS:
%   Z := Point of evaluation (p,u) given as a numeric matrix (IE [6 x 3]
%       for TSO(3)
%   X, Y := Vector fields on TSO(n) given as numeric matrices (IE: [6 x 3]
%       for TSO(3)
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   D := Result of <X,Y> on TSO(n) 
%% Standard Variables
if (length(m) ~= 3)
    error('m must be an array of length 3')
end
m1 = m(1); m2 = m(2); m3 = m(3); % These should produce a pos. def. symm. matrix M

%% Split X and Y into vertical and horizontal components
[R, ~] = rotBundle_split(Z);
Xh = rotBundle_extractHoriz(Z, X);
Xv = rotBundle_extractVert(Z, X);
Yh = rotBundle_extractHoriz(Z, Y);
Yv = rotBundle_extractVert(Z, Y);

%% <Xh, Yh> on TSO(n) == m1<X,Y> on SO(n)
metric_Xh_Yh = m1*rot_metric(R, Xh, Yh);

%% <Xv, Yh> on TSO(n) == m2<X,Y> on SO(n)
metric_Xv_Yh = m2*rot_metric(R, Xv, Yh);
metric_Xh_Yv = m2*rot_metric(R, Xh, Yv);

%% <Xv, Yv> on TSO(n) == m3<X,Y> on SO(n)
metric_Xv_Yv = m3*rot_metric(R, Xv, Yv);

%% Resulting D is the sum of the three components
D = metric_Xh_Yh + metric_Xv_Yh + metric_Xh_Yv + metric_Xv_Yv;
end

