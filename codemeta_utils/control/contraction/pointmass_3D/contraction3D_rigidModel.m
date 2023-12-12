function [ dx ] = contraction3D_rigidModel( x, m, F )
%Create a rigid point mass model in R^3. Model is given by m*ddx=F
% INPUT:
% x - the starte vector [3x1]
% v - the velocity vector [3x1]
% m - mass
% F - force

dx = [zeros(3,3) eye(3);zeros(3,3) zeros(3,3)]*x + [zeros(3,3);eye(3)]*F*1/m;
end

