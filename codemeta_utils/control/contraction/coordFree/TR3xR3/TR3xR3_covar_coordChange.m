function [D_X_Y_nonnatural, D_X_Y_natural] = TR3xR3_covar_coordChange(A,B,C,t,X,Y,M_nn,varargin)
% Compute the Levi-Civita connection coompatible with the nonnatural metric
% M_nn on TR^3xR^3 at (A,B,C)|_t=t by change of coordinates from nonnatural
% to natural metric
% INPUTS:
%   A(t) := The evaluation position on TR^3
%   B(A,t) := The evaluation velocity on TR^3
%   C(t) := The evaluation position on R^3
%   t := Time of evaluation
%   X(A,B,C), Y(A,B,C) := Tangent Vector as [9x1] vector
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.
% OOUTPUTS:
%   D_X_Y_nonnatural, D_X_Y_natural := The covariant derivative of Y along 
%       X evaluated at (A,B,C)|_t=t.

% Optional Parameters
flagUseNewCurve = false;
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'usenewcurve'
            % Flag to use the new curve resulting from the transformed
            % vector field
            flagUseNewCurve = true;
        otherwise    
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

% Evaluate the point
A_t = A(t); B_t = B(A(t),t); C_t = C(t);
[J,M_n] = TR3xR3_SchurComplement(M_nn);
% Extract the 3 tangent vectors of each vector field
X_pos = @(A,B,C) extractComp(X(A,B,C),1,3,1,1);
X_vel = @(A,B,C) extractComp(X(A,B,C),4,6,1,1);
X_ref = @(A,B,C) extractComp(X(A,B,C),7,9,1,1);

% Perform coordinate transform
% X_natural = @(A,B,C) kron(J,eye(3))*X(A,B,C);
Y_natural = @(A,B,C) kron(J,eye(3))*Y(A,B,C);


% Use X_natural to define new curve to diff along (THIS IS WHAT WE'RE DOING
% ON TSO3XSO3, BUT RESULT SHOULD BE TO USE ORIGINAL CURVE AND NOT THE NEW
% ONE)
if flagUseNewCurve
    A_curve = @(t) A_t + X_pos(A_t,B_t,C_t)*t;
    B_curve = @(A,t) B_t + X_vel(A_t,B_t,C_t)*t;
    C_curve = @(t) C_t + (X_ref(A_t,B_t,C_t) + J(3,1)*X_pos(A_t,B_t,C_t) + J(3,2)*X_vel(A_t,B_t,C_t))*t;
else
    A_curve = A;
    B_curve = B;
    C_curve = C;
end

D_X_Y_natural = funApproxDer(@(t) Y_natural(A_curve(t),B_curve(A_curve(t),t),C_curve(t)),t);
D_X_Y_nonnatural = kron(inv(J),eye(3))*D_X_Y_natural;
end



