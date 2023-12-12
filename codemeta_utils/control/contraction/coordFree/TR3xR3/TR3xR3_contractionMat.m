function [M_contract] = TR3xR3_contractionMat(kd,kv,M_nn)
% Compute the matrix form of the M matrix (on TR^3xR^3) such that
% [zeta;eta;nu]'*M*[zeta;eta;nu] = <D_Y_X+lambda*Y,Y> where zeta and eta
% are the R^3 representation of the horizontal and vertical tangent vectors
% of Y on TR^3 and nu is the R^3 vector on R^3. Assume cost function is
% f(x,xd) = 1/2*norm(x-xd)^2
% INPUTS:
%   kd,kv := scalar position/velocity gains
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.
% OUTPUTS:
%   M_contract := M matrix such that 
%       [zeta;eta;nu]'*M*[zeta;eta;nu] = <D_Y_X+lambda*Y,Y>
%       NOTE: This should be a function of (A,B,C) but it doesnt matter in
%       R^6

% Extract the gains
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2);
m4 = M_nn(3,3); m5 = M_nn(2,3); m6 = M_nn(1,3);

% Define the components
M11 = -kd*m2*eye(3);
M21 = eye(3)*m1/2 - m3*kd*eye(3)/2 - kv*m2*eye(3)/2;
M22 = eye(3)*m2 - kv*m3*eye(3);
M31 = -eye(3)*m6/2 - m5*kd*eye(3)/2 + kd*m2*eye(3)/2;
M32 = eye(3)*m6/2 - eye(3)*m5/2 - kv*m5*eye(3)/2 + kd*m3*eye(3)/2;
M33 = -eye(3)*m4 + kd*m5*eye(3);

M_contract = [M11 M21' M31';M21 M22 M32';M31 M32 M33];
end

