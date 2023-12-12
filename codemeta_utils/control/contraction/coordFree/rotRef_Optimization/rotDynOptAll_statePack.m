function [x] = rotDynOptAll_statePack(R,w,RReference,kd,kv,kref,M_nn,b)
% Last Edited: Nov 12, 2020 by Bee Vang
% Pack the states given as matrices into a column vector for the numerical
% integrator
% INPUTS:
%   R := [3x3] rotation matrix of current position
%   w := [3x1] current velocity vector
%   RReference := [3x3] rotation matrix of current refence position
%   kd := positive scalar proportional gain
%   kv := positive scalar deriviative gain
%   kref := positive scalar
%   M_nn := [3x3] non-natural contraction metric
%   b := positive scalar for min. exp. conv. rate bound
% OUTPUTS:
%   x := [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]

x(1:9,1) = R(:);
x(10:12,1) = w(:);
x(13:21,1) = RReference(:);
x(22,1) = kd;
x(23,1) = kv;
x(24,1) = kref;
x(25:33,1) = M_nn(:);
x(34,1) = b;
end

