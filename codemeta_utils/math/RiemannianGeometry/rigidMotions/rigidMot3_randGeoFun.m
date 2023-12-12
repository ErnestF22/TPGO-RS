function [Gt,dGt,G0,dG0,vGVec,ddGt,dvGVec] = rigidMot3_randGeoFun(A,varargin)
% Randomly generate a geodesic curve on SE(3)
% If A is given, generate a geodesic starting at A
% INPUTS:
%   A := a position on SE(3)
% OUTPUTS:
%   Gt := Geodesic curve as a function of time
%   dGt := Velocity along the curve Gt
%   G0 := Initial position on SE(3)
%   dG0 := Initial tangent vector
%   vVec := The R^6 vector representation of dGt
%   ddGt := Jerk along the curve Gt
%   dvVec := The R^6 vector representation of ddGt

% Parameters
v = cnormalize(randn(3,1)); % Velocity vector of geodesic on R^3

if ~exist('A', 'var') || isempty(A)
    % Generate a random rotation and position
    R0 = rot_randn; d0 = randn(3,1);
    A = [R0, d0;0 0 0 1];
end

%  The shortest path on SE(3) is obtained by lifiting the geodesics from
%  SO(3) and R^3.
R0 = A(1:3,1:3); d0 = A(1:3,4);
[Rt, dRt, ~, dR0, vVec, ddRt, dvVec]  = rot_randGeodFun(R0);
dt = @(t) t*v+d0;
vt = @(t) v;
% Combine the geodesics from SO(3) and R^3
Gt = @(t) [Rt(t) dt(t);0 0 0 1];
dGt = @(t) [dRt(t) vt(t);0 0 0 1];
G0 = A;
dG0 = [dR0 v;0 0 0 1];
vGVec = [vVec;v];
ddGt = @(t) [ddRt(t) zeros(3,1);0 0 0 1];
dvGVec = [dvVec(0); zeros(3,1)];
end

