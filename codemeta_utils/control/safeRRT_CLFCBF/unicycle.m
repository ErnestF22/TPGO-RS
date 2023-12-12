function [dx] = unicycle(x,u)
% First order dynamical model of a unicycle
% INPUTS:
%   x := current state vector [3x1]
%   u := control vector [2x1]
% OUTPUTS:
%   dx := updated dynamics

x1 = x(1); x2 = x(2); theta = x(3);
v = u(1); w = u(2);

dx(1,1) = cos(theta)*v;
dx(2,1) = sin(theta)*v;
dx(3,1) = w;

end

