% Test vector fields on TSO3 with a componenet dependent on a reference SO3
% manifold (by left translation)

close all; clear all; clc;
fprintf('Test Script: "%s"\n\n',mfilename)
%% Define Parameters
t = rand; % Pick a random time to evaluate
% Define curves on TSO(3)
[R,~,R0,dR0] = rot_randGeodFun;
dR = @(R) R*R0'*dR0;
% Linear U
uVec = randn(3,1);
U = @(R,t) R*hat3(t*uVec);
du = uVec;
% Define curve on reference manifold SO(3)
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Random non-natural metric on TSO(3)
M_nonnatural = randn(2,2); M_nonnatural = M_nonnatural'*M_nonnatural;
m_nonnatural = [M_nonnatural(1,1);M_nonnatural(1,2);M_nonnatural(2,2)];
% Define Left Invarient Vector Fields on TSO(3)
x1 = randn(3,1); x2 = randn(3,1);
X = @(R,U,RRef) [U;R*hat3(x1)+R*RRef'*RRef*hat3(x2)];
% X = @(R,U,RRef) [U;R*hat3(x1)+R*hat3(x2)];
X_R = @(R,U) X(R,U,RRef(t));
% Generate vector field along the curve Z(t) = [R(t);U(R,t)]
Z = @(R,U) [R;U];
Y = @(R,U) [dR(R);R*hat3(du)]; % NOTE: dZ/dt = Y;
U_R = @(R) U(R,t);
LT_R = @(R) kron(eye(2),R); % Convient matrix for Left translation

%% Compute the change of the nonnatural metric (ie d/dt<X,Y>)
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t));
Zt = @(t) Z(R(t),U(R(t),t));

g_nonnatural = @(t) rotBundle_metric_nonNatural(Zt(t),Xt(t),Yt(t),m_nonnatural);
g_nonnatural_der = funApproxDer(g_nonnatural,t);
%% Compute the change of the nonnatural metric using covar ders (ie <D_Y_X,Y> + <D_Y_Y,X>)

% Compute D_Y_X
D_Y_X = rotBundle_covar_nonNatural(Y,X_R,U_R,m_nonnatural,'rigidrot',0);
% Compute <D_Y_X,Y>
g_DYX_Y = rotBundle_metric_nonNatural(Zt(t),D_Y_X(R(t)),Yt(t),m_nonnatural);

% Compute D_X_Y
D_Y_Y = rotBundle_covar_nonNatural(Y,Y,U_R,m_nonnatural);
% Compute <D_Y_Y,X>
g_DYY_X = rotBundle_metric_nonNatural(Zt(t),D_Y_Y(R(t)),Xt(t),m_nonnatural);

g_analytical_der = g_DYX_Y + g_DYY_X;
%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('analytical der: %0.5f:\n',g_analytical_der);
fprintf('Error of der: %0.5f:\n',g_nonnatural_der-g_analytical_der);