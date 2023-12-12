function [R,w,RReference,kd,kv,kref,M_nn,b] = rotDynOptAll_stateUnpack(x)
% Last Edited: Nov 12, 2020 by Bee Vang
% Extract the individual elements of the state vector
% x = [R,w,RRef,kd,kv,kref,M_nn,beta]
% INPUTS:
%   x := [34x1] vector with states [R,w,RRef,kd,kv,kref,M_nn,beta]
% OUTPUTS:
%   R := [3x3] rotation matrix of current position
%   w := [3x1] current velocity vector
%   RReference := [3x3] rotation matrix of current refence position
%   kd := positive scalar proportional gain
%   kv := positive scalar deriviative gain
%   kref := positive scalar
%   M_nn := [3x3] non-natural contraction metric
%   b := positive scalar for min. exp. conv. rate bound

% Extract states
[R,w,RReference] = rotDyn_stateUnpack(x,'augmentedsystem');
% Extract gains
kd = x(22,:);
kv = x(23,:);
kref = x(24,:);
% Extract metric
M_nn = reshape(x(25:33,:),3,3,[]);
% Extract convergence rate
b = x(34,:);
end