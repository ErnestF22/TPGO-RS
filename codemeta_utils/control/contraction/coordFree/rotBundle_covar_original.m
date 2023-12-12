function [ nablaY_hat ] = rotBundle_covar_original( X, Y)
% Generate a function handle for the connection on the tangent bundle of
% two vector fields on TSO(n)
% INPUTS:
%   X(R,U), Y(R,U) := vector fields on TSO(n) represented by function 
%       handles so that X(R,U) returns the vector field X/Y evaluated 
%       at (R,U).
% OUTPUTS:
%   nablaY_hat(R,U) := the connection on TSO(n) of Y with respect to X
%       represented as a function of R

% The connection on TSO(n) is the sum of four components given by Prop. 7.2
% of Gudmundsson and Kappos (pg. 15). 

%% General use objects
% The point (p,u) of interest on the tangent bundle TSO(n)
Z = @(R, U) [R;U];
% The horizontal components of X and Y as tangent vectors in T_{p}SO(n)
% Xh = @(R,U) rotBundle_extractHoriz(Z(R,U), X(R,U));
Yh = @(R,U) rotBundle_extractHoriz(Z(R,U), Y(R,U));
% The vertical components of X and Y as tangent vectors in T_{p}SO(n)
Xv = @(R,U) rotBundle_extractVert(Z(R,U), X(R,U));
% Yv = @(R,U) rotBundle_extractVert(Z(R,U), Y(R,U));

% Connection of Y^{h} wrt X^{h} on TSO(n)
DelHat_Xh_Yh = @(R,U) funDelHat_Xh_Yh(R, U, X, Y);
% Connection of Y^{v} wrt X^{h} on TSO(n)
DelHat_Xh_Yv = @(R,U) funDelHat_Xh_Yv(R, U, X, Y);
% Connection of Y^{h} wrt X^{v} on TSO(n)
DelHat_Xv_Yh = @(R,U) funDelHat_Xv_Yh(R, U, X, Y);
% Connection of Y^{v} wrt X^{v} on TSO(n)
DelHat_Xv_Yv = @(R,U) funDelHat_Xv_Yv(R, U, X, Y);
% Return the connection of Y wrt X at Z on TSO(n)
nablaY_hat = @(R,U) DelHat_Xh_Yh(R,U) + DelHat_Xh_Yv(R,U) + ...
    DelHat_Xv_Yh(R,U) + DelHat_Xv_Yv(R,U);
end

function [ hh ] = funDelHat_Xh_Yh(R, U, X, Y)
% Given R, U find the connection of Y^{h} wrt X^{h} on TSO(n)
% NOTE: everything in this function is happening on SO(n) thus components 
% of X, Y are functions of R only.
% INPUTS:
%   R, U := The point of interest on TSO(n) as numerical matrices.
%   X(R,U), Y(R,U) := Tangent vectors \in T_{R}T_{U}SO(n) as function
%       handles of paramters (R,U)
% OUTPUTS:
%   hh := The numerical tangent vector \in T_{R}T_{U}SO(3) of the
%       connection of Y^{h} wrt X^{h}.

% Define the point of interest on TSO(n)
Z = [R; U];
% Extract the horizontal components of the X, Y tangent vectors as tangent 
% vectors in T_{R}SO(n)
Xh = rotBundle_extractHoriz(Z, X(R,U));
Yh = @(R) rotBundle_extractHoriz(Z, Y(R,U));
% Find the connection of Yh wrt Xh on SO(3)
gammaX = @(t) rot_exp(R, t*Xh); % A path that produces Xh
dY = funApproxDer(@(t) Yh(gammaX(t)), 0); % Time derivative of Yh along gammaX
Del_Xh_Yh = rot_covar(Xh, Yh(R), dY, R);
% Horizontal lift of the resulting connection
Del_Xh_Yh_h = rotBundle_horizLift(Del_Xh_Yh);
% Find the curvature tensor (R(X,Y)U)
Rp = rot_curvature(R, Xh, Yh(R), U);
% Vertical lift of Rp
Rp_v = rotBundle_vertLift(Rp);
%
hh = Del_Xh_Yh_h - 1/2*Rp_v;
end

function [ hv ] = funDelHat_Xh_Yv(R, U, X, Y)
% Given R, U find the connection of Y^{v} wrt X^{h} on TSO(n)
% NOTE: everything in this function is happening on SO(n) thus components 
% of X, Y are functions of R only.
% INPUTS:
%   R, U := The point of interest on TSO(n) as numerical matrices.
%   X(R,U), Y(R,U) := Tangent vectors \in T_{R}T_{U}SO(n) as function
%       handles of paramters (R,U)
% OUTPUTS:
%   hv := The numerical tangent vector \in T_{R}T_{U}SO(3) of the
%       connection of Y^{v} wrt X^{h}.

% Define the point of interest on TSO(n)
Z = @(R) [R; U];
% Extract the horizontal component of the X vector as tangent a
% tangent vector in T_{R}SO(n)
Xh = rotBundle_extractHoriz(Z(R), X(R,U));
% Extract the vertical component of the Y vector as tangent a
% tangent vector in T_{R}SO(n)
Yv = @(R) rotBundle_extractVert(Z(R), Y(R,U));
% Find the connection of Yv wrt Xh on SO(3)
gammaX = @(t) rot_exp(R, t*Xh); % A path that produces Xh
dY = funApproxDer(@(t) Yv(gammaX(t)), 0); % Time derivative of Yv along gammaX
Del_Xh_Yv = rot_covar(Xh, Yv(R), dY, R);
% Horizontal lift of the resulting connection
Del_Xh_Yh_v = rotBundle_vertLift(Del_Xh_Yv);
% Find the curvature tensor (R(U,Y)X)
Rp = rot_curvature(R, U, Yv(R), Xh);
% Horizontal lift of Rp
Rp_h = rotBundle_horizLift(Rp);
%
hv = Del_Xh_Yh_v + 1/2*Rp_h;
end

function [ vh ] = funDelHat_Xv_Yh(R, U, X, Y)
% Given R, U find the connection of Y^{h} wrt X^{v} on TSO(n)
% NOTE: everything in this function is happening on SO(n) thus components 
% of X, Y are functions of R only.
% INPUTS:
%   R, U := The point of interest on TSO(n) as numerical matrices.
%   X(R,U), Y(R,U) := Tangent vectors \in T_{R}T_{U}SO(n) as function
%       handles of paramters (R,U)
% OUTPUTS:
%   vh := The numerical tangent vector \in T_{R}T_{U}SO(3) of the
%       connection of Y^{h} wrt X^{v}.

% Define the point of interest on TSO(n)
Z = [R;U];
% Extract the vertical component of the X vector as tangent a
% tangent vector in T_{R}SO(n)
Xv = rotBundle_extractVert(Z, X(R,U));
% Extract the horizontal component of the Y vector as tangent a
% tangent vector in T_{R}SO(n)
Yh = rotBundle_extractHoriz(Z, Y(R,U));
% curvature tensor (R(u,X)Y)
Rp = rot_curvature(R, U, Xv, Yh);
% Horizontal lift of Rp
Rp_h = rotBundle_horizLift(Rp);
% The connection of Y^{h} wrt X^{v} == 1/2*Rp_h
vh = 1/2*Rp_h;
end

function [ vv ] = funDelHat_Xv_Yv(R, U, X, Y)
% Given R, U find the connection of Y^{v} wrt X^{v} on TSO(n)
% NOTE: everything in this function is happening on SO(n) thus components 
% of X, Y are functions of R only.
% INPUTS:
%   R, U := The point of interest on TSO(n) as numerical matrices.
%   X(R,U), Y(R,U) := Tangent vectors \in T_{R}T_{U}SO(n) as function
%       handles of paramters (R,U)
% OUTPUTS:
%   vh := The numerical tangent vector \in T_{R}T_{U}SO(3) of the
%       connection of Y^{v} wrt X^{v}.

% The connection of Y^{v} wrt X^{v} is zero
vv = zeros(size([R;U]));
end