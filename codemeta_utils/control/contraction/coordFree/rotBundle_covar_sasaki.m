 function [ nablaY_hat ] = rotBundle_covar_sasaki( X, Y, U)
% Generate a function handle for the connection on the tangent bundle of
% two vector fields on TSO(n)
% INPUTS:
%   X(R,U), Y(R,U) := vector fields on TSO(n) represented by function 
%       handles so that X(R,U) returns the vector field X/Y evaluated 
%       at (R,U).
%   U(R) := the tangent vector at R which creates the curve on TSO(n) in
%       the form of [R; U]
% OUTPUTS:
%   nablaY_hat(R) := the connection on TSO(n) of Y with respect to X
%       represented as a function of R

% The connection on TSO(n) is the sum of four components given by Prop. 7.2
% of Gudmundsson and Kappos (pg. 15). 

%% General use objects
% The point (p,u) of interest on the tangent bundle TSO(n)
Z = @(R) [R;U(R)];
% The horizontal components of X and Y as tangent vectors in T_{p}SO(n)
Xh = @(R) rotBundle_extractHoriz(Z(R), X(R,U(R)));
Yh = @(R) rotBundle_extractHoriz(Z(R), Y(R,U(R)));
% The vertical components of X and Y as tangent vectors in T_{p}SO(n)
Xv = @(R) rotBundle_extractVert(Z(R), X(R,U(R)));
Yv = @(R) rotBundle_extractVert(Z(R), Y(R,U(R)));

%% Connection of Y^{h} wrt X^{h} on TSO(n)
% Connection of Yh wrt Xh on the manifold
% gammaX = @(R,t) rot_exp(R,t*Xh(R));
% dY1 = @(R) funApproxDer( @(t) Yh(gammaX(R,t)), 0);
% Del_Xh_Yh = @(R) rot_covar(Xh(R), Yh(R), dY1(R), R);
Del_Xh_Yh = rot_covar(Xh, Yh);
% horizontal lift of Del_X_Y(R)
Del_X_Y_h = @(R) rotBundle_horizLift(Del_Xh_Yh(R));
% curvature tensor (R(X,Y)U)
Rp1 = @(R) rot_curvature(R, Xh(R), Yh(R), U(R));
% vertical lift of Rp1
Rp1_v = @(R) rotBundle_vertLift(Rp1(R));

DelHat_Xh_Yh = @(R) Del_X_Y_h(R) - 1/2*Rp1_v(R);
%% Connection of Y^{v} wrt X^{h} on TSO(n)
% Connection of Yv wrt Xh on the manifold
% dY2 = @(R) funApproxDer( @(t) Yv(gammaX(R,t)), 0);
% Del_Xh_Yv = @(R) rot_covar(Xh(R), Yv(R), dY2(R), R);
Del_Xh_Yv = rot_covar(Xh, Yv);
% vertical lift of Del_X_Y(R)
Del_X_Y_v = @(R) rotBundle_vertLift(Del_Xh_Yv(R));
% curvature tensor (R(u,Y)X)
Rp2 = @(R) rot_curvature(R, U(R), Yv(R), Xh(R));
% horizontal lift of Rp2
Rp2_h = @(R) rotBundle_horizLift(Rp2(R));

DelHat_Xh_Yv = @(R) Del_X_Y_v(R) + 1/2*Rp2_h(R);
%% Connection of Y^{h} wrt X^{v} on TSO(n)
% curvature tensor (R(u,X)Y)
Rp3 = @(R) rot_curvature(R, U(R), Xv(R), Yh(R));
% horizontal lift of Rp3
Rp3_h = @(R) rotBundle_horizLift(Rp3(R));

DelHat_Xv_Yh = @(R) 1/2*Rp3_h(R);
%% Connection of Y^{v} wrt X^{v} on TSO(n)
DelHat_Xv_Yv = @(R) zeros(6,3);

%% Return the connection of Y wrt X at Z on TSO(n)
nablaY_hat = @(R) DelHat_Xh_Yh(R) + DelHat_Xh_Yv(R) + ...
    DelHat_Xv_Yh(R) + DelHat_Xv_Yv(R);
end

