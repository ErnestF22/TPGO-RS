function [CLF,Grad_CLF] = compute_CLF(x,xd)
% Compute the CLF and its partial derivatives
% INPUTS:
%   x := current state vector [3x1]
%   xd := desired state [3x1]
% OUTPUTS:
%   CLF := the CLF evaluated at (x,xd)
%   Grad_CLF := partial derivative of CLF wrt [x1;x2;x3], this is a [1x3]
%       vector

% Extract/Redefine variables
x1 = x(1); x2 = x(2); x3 = x(3);
x1s = xd(1); x2s = xd(2); x3s = xd(3);

CLF = (x1 - x1s)^2 + (x2 - x2s)^2 + (cos(x3 - x3s) - (x1 - x1s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))^2 + (sin(x3 - x3s) - (x2 - x2s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))^2;
Grad_CLF = [2*x1 - 2*x1s - 2*(cos(x3 - x3s) - (x1 - x1s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))*(1/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2) - (abs(x1 - x1s)*sign(x1 - x1s)*(x1 - x1s))/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(3/2)) + (2*abs(x1 - x1s)*sign(x1 - x1s)*(sin(x3 - x3s) - (x2 - x2s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))*(x2 - x2s))/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(3/2),...
    2*x2 - 2*x2s - 2*(sin(x3 - x3s) - (x2 - x2s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))*(1/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2) - (abs(x2 - x2s)*sign(x2 - x2s)*(x2 - x2s))/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(3/2)) + (2*abs(x2 - x2s)*sign(x2 - x2s)*(cos(x3 - x3s) - (x1 - x1s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))*(x1 - x1s))/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(3/2),...
    2*cos(x3 - x3s)*(sin(x3 - x3s) - (x2 - x2s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2)) - 2*sin(x3 - x3s)*(cos(x3 - x3s) - (x1 - x1s)/(abs(x1 - x1s)^2 + abs(x2 - x2s)^2)^(1/2))];
end

