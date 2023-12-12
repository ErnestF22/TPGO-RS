function [V,Vdot] = rotBundle_LyapunovW(w, kv)
% Compute the Lyapunov function for the velocity only controller
% INPUTS:
%   w := Current velocity (array of [3x1] vectors)
%   kv := Velocity error gain (scalar)
% OUTPUTS:
%   V := Lyapunov value(s) (should be >= 0)
%   Vdot := Derivative of V

V = vecnorm(w);
Vdot = -kv;
% V = vecnorm(w).^2;
% Vdot = -kv*V;

end