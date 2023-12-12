% Compute and compare the contraction metric when computed using the direct
% method and coordinate change method.
close all; clear all; clc;

%% Define Parameters
syms m1 m2 m3 m4 m5 m6 kd kv real
zeta = sym('zeta',[3 1],'real');
eta = sym('eta',[3 1],'real');
nu = sym('nu',[3 1],'real');
Y = [zeta;eta;nu];
M_nn = [m1 m2 m6;m2 m3 m5;m6 m5 m4];
[J,M_n]=TR3xR3_SchurComplement(M_nn);

%% Direct method
D_Y_X_direct = [eta;-kd*zeta-kv*eta+kd*nu;-nu];
g_DYX_Y_direct = D_Y_X_direct'*kron(M_nn,eye(3))*Y;
g_DYX_Y_direct_analytical = eta'*m1*zeta + eta'*m2*eta + eta'*m6*nu...
    +(-kd*zeta'-kv*eta'+kd*nu')*(m2*zeta + m3*eta + m5*nu)...
    -nu'*m6*zeta-nu'*m5*eta-nu'*m4*nu;
fprintf('g_DYX_Y_direct-g_DYX_Y_direct_analytical: %s\n',...
    char(simplify(g_DYX_Y_direct-g_DYX_Y_direct_analytical)));
%% Coord change method (using new JY as vector field to diff along)
A = -kd*zeta-kv*eta+kd*(nu+m6/m4*zeta+m5/m4*eta);
B = -(nu+m6/m4*zeta+m5/m4*eta) + m6/m4*eta...
    -m5/m4*kd*zeta + m5/m4*kd*(nu+m6/m4*zeta+m5/m4*eta)...
    -m5/m4*kv*eta;
D_Y_X_coord = [eta;A;B];

g_DYX_Y_coord = D_Y_X_coord'*kron(M_n,eye(3))*(kron(J,eye(3))*Y);

%% Coord change method (using original "nonnatural" vector field to dif along)
D_Y_X_coord_2 = [eta;...
    -kd*zeta-kv*eta+kd*nu;...
    -nu+m6/m4*eta-m5/m4*kd*zeta+m5/m4*kd*nu-m5/m4*kv*eta];
g_DYX_Y_coord_2 = D_Y_X_coord_2'*kron(M_n,eye(3))*(kron(J,eye(3))*Y);

%% Show results
fprintf('Metric Error when using new JY to define curve:\n');
simplify(g_DYX_Y_direct-g_DYX_Y_coord)
fprintf('Metric Error when using original nonnatural curve:\n');
simplify(g_DYX_Y_direct-g_DYX_Y_coord_2)
