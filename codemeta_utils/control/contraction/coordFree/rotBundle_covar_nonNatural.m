function [ nablaY_hat ] = rotBundle_covar_nonNatural( X, Y, U, m, varargin)
% Mult <D_X_Y, M*Z> by M^{-1}, then set Zh = 0 to get vertical component,
% and set Zv = 0 to get horizontal component

% Compute the connection on the tangent bundle of two vector fields on TSO(3)
% at (R,U(R))|_{t=t}. 
% NOTE: Currently this is returned as a function of R, but R must be
% evaluated at R|_{t=t} otherwise the result does not hold
% INPUTS:
%   X(R,U), Y(R,U) := vector fields on TSO(n) represented by function 
%       handles so that X(R,U) returns the vector field X/Y evaluated 
%       at (R,U).
%   U(R) := the tangent vector at R which creates the curve on TSO(n) in
%       the form of [R; U]
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   nablaY_hat(R) := the connection on TSO(n) of Y with respect to X
%       represented as a function of R

% The connection on TSO(n) is the sum of four components given by Prop. 7.2
% of Gudmundsson and Kappos (pg. 15). 

%% General use objects
if (length(m) ~= 3)
    error('m must be an array of length 3')
end
m1 = m(1); m2 = m(2); m3 = m(3); % These should produce a pos. def. symm. matrix M
flagRigidRotControl = false; % Assumes Y is the controled closed loop VF for rigid body rotations
kv = 0;
%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'rigidrot'
            % Add adjustment terms since Y has some dependency on its
            % position along the fiber. IE Y is not completely
            % lift-decomposable.
            % Assumes VF is rigid body rotation with controller 
            % u = kd*log(R,eye(3)) - kv*Yh
            flagRigidRotControl=true;
            % Must also send in a value for kv gain!
            ivarargin=ivarargin+1;
            kv = varargin{ivarargin};
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 < len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end
    end
    ivarargin=ivarargin+1;
end

% The point (p,u) of interest on the tangent bundle TSO(n)
Z = @(R) [R;U(R)];
% The horizontal components of X and Y as tangent vectors in T_{p}SO(n)
Xh = @(R) rotBundle_extractHoriz(Z(R), X(R,U(R)));
Yh = @(R) rotBundle_extractHoriz(Z(R), Y(R,U(R)));
% The vertical components of X and Y as tangent vectors in T_{p}SO(n)
Xv = @(R) rotBundle_extractVert(Z(R), X(R,U(R)));
Yv = @(R) rotBundle_extractVert(Z(R), Y(R,U(R)));

% Connection of Yh wrt Xh on the manifold
Del_Xh_Yh = rot_covar(Xh, Yh);
% curvatur tensor (R(U,X)Y)
R_U_Xh_Yh = @(R) rot_curvature(R, U(R), Xh(R), Yh(R));
% curvature tensor (R(X,Y)U)
R_Xh_Yh_U = @(R) rot_curvature(R, Xh(R), Yh(R), U(R));
% Connection of Yv wrt Xh on the manifold
Del_Xh_Yv = rot_covar(Xh, Yv);
% curvature tensor (R(U,Y)X)
R_U_Yv_Xh = @(R) rot_curvature(R, U(R), Yv(R), Xh(R));
% curvature tensor (R(U,X)Y)
R_U_Xv_Yh = @(R) rot_curvature(R, U(R), Xv(R), Yh(R));

%% Define the horizontal components
Del_bar_X_Y_h = @(R) rotBundle_horizLift(...
    (m1*m3*Del_Xh_Yh(R) + m2*m3*R_U_Xh_Yh(R)...
    - m2^2*Del_Xh_Yh(R) + m2*m3/2*R_Xh_Yh_U(R)...
    + m3^2/2*R_U_Yv_Xh(R)...    
    + m3^2/2*R_U_Xv_Yh(R))/(m1*m3-m2^2));
%% Define the vertical components
Del_bar_X_Y_v = @(R) rotBundle_vertLift(...
    (- m2^2*R_U_Xh_Yh(R)...
    - m1*m3/2*R_Xh_Yh_U(R)...
    - m2^2*Del_Xh_Yv(R) - m2*m3/2*R_U_Yv_Xh(R)...
    + m1*m3*Del_Xh_Yv(R)...
    - m2*m3/2*R_U_Xv_Yh(R))/(m1*m3-m2^2));

%% Return the connection of Y wrt X at Z on TSO(n)
nablaY_hat = @(R) Del_bar_X_Y_h(R) + Del_bar_X_Y_v(R);

if flagRigidRotControl
    % add terms to account for changing vector field along fibers
    nablaY_hat = @(R) nablaY_hat(R) + [Xv(R);-kv*Xv(R)];
end

