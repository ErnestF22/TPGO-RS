function [D_X_Y_nonnatural] = TR3xR3_covar_directMethod(A,B,C,t,X,Y,M_nn)
% Compute the Levi-Civita connection coompatible with the nonnatural metric
% M_nn on TR^3xR^3 at (A,B,C)|_t=t by direct method (knowning how to
% compute covar on R^6)
% INPUTS:
%   A(t) := The evaluation position on TR^3
%   B(A,t) := The evaluation velocity on TR^3
%   C(t) := The evaluation position on R^3
%   t := Time of evaluation
%   X(A,B,C), Y(A,B,C) := Tangent Vector as [9x1] vector
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.
% OOUTPUTS:
%   D_X_Y_nonnatural := The covariant derivative of Y along X evaluated at
%   (A,B,C)|_t=t.

% % Evaluate the point
% A_t = A(t); B_t = B(t); C_t = C(t);
% % Redefine B as function of A
% B_A = @(A) B(A,t);

% NOTE: to be more correct we should take the derivative of Y along any
% curve with tangent vector X... (we assume that X = d/dt along [A;B;C])
D_X_Y_nonnatural = funApproxDer(@(t) Y(A(t),B(A(t),t),C(t)),t);
end

