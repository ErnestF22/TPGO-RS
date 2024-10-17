function [nablaY, nablaYVec] = rigidMot3_covar(X,Y)
% Compute the covariant derivative of the vector field Y along vector field
% X on SE(3).
% INPUTS:
%   X,Y := Vector fields on SE(3) as functions of any curve G=fun(R,d). 
%       The VF are represented as [4x4] matrices.
% OUTPUTS:
%   nablaY := The covarient derivative of the vector field Y with respect
%       to vector field X as a function of curve G
%       nablaYVec := Vectorize version of nablaY

if ~isa(Y, 'function_handle')
    error('Numerics are not implemented yet');
end

% The connection on SE(3) using the standard metric eqn (7) is given in
% "Planning of smooth motions on SE(3)" (see Zefran and Kumar) in eqn (10).

% Extract the rotation and position related tangent vectors from X,Y
X_R = @(R,d) rigidMot3_extractRot(X(R,d));
X_RVec = @(R,d) rot_vee(R,X_R(R,d));
X_d = @(R,d) rigidMot3_extractPos(X(R,d));
Y_R = @(R,d) rigidMot3_extractRot(Y(R,d));
% Y_RVec = @(R,d) rot_vee(R,X_R(R,d));
Y_d = @(R,d) rigidMot3_extractPos(Y(R,d));

% Take the numerical time derivative of Y along any curve with tangent
% vector X (for simplicity, a geodesic)
gammaR_X = @(R,d,t) rot_exp(R, t*X_R(R,d)); % great circle on SO(3)
gammaR_d = @(R,d,t) t*X_d(R,d) + d; % line in R^3
dY_R = @(R,d) funApproxDer( @(t) Y_R(gammaR_X(R,d,t),gammaR_d(R,d,t)),0); %This is not a tangent vector,
                % how to vectorize as dw_{y}/dt (seems that this is rot_covar... wtih some dependence on position)
dY_d = @(R,d) funApproxDer( @(t) Y_d(gammaR_X(R,d,t),gammaR_d(R,d,t)),0); %This is dv_{y}/dt

% Construct the covar according to Zefran and Kumar
% Use Edelman results for covar on SO(3)
Gamma_R = @(R,d) R*(X_R(R,d)'*Y_R(R,d) + Y_R(R,d)'*X_R(R,d))/2;
nablaY_R = @(R,d) dY_R(R,d) + Gamma_R(R,d);
nablaY_RVec = @(R,d) rot_vee(R,nablaY_R(R,d));
% Use Zefran & Kumar
nablaY_d = @(R,d) dY_d(R,d) + cross(X_RVec(R,d), Y_d(R,d));

% Results
nablaY = @(R,d) [nablaY_R(R,d) nablaY_d(R,d);0 0 0 1];
nablaYVec = @(R,d) [nablaY_RVec(R,d); nablaY_d(R,d)];
end

