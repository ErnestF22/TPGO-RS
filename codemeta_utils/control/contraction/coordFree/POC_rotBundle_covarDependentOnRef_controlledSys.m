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
% Define a closed looped PD controlled attitude system on TSO3 (essentially
% treat the RRef dependent component as left-invar till we take the total
% time derivative by chain rule)
X = @(R,U,RRef) [U+.05*R*RRef'*rot_log(RRef,eye(3));5*rot_log(R,RRef)-2*U+R*RRef'*rot_log(RRef,eye(3))];
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
% Define useful constants
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t); % Evaluated at t

% Compute D_Y_X
D_Y_X = rotBundle_covar_nonNatural(Y,X_R,U_R,m_nonnatural,'rigidrot',2);
% Additional term to account for the tangent vector from RRef depending on
% RRef (this is seen as a direct time dependence on TSO3)
dX_dt = funApproxDer(@(t) X(R_t,U_t,RRef(t)),t);
D_Y_X = @(R) D_Y_X(R) + dX_dt;
M_test = zeros(3,3); M_test(2:3,2:3) = M_nonnatural;
D_Y_X2 = rotRef_covar(@(R,U,RRef) [zeros(3);Y(R,U)],....
    @(R,U,RRef) [zeros(3); X(R,U,RRef)],...
    R,U_R,RRef,M_test,t,'rigidrot',2);
fprintf('Test rotRef_covar (should be zero), D_Y_X at (R,U):\n');
D_Y_X2(R(t)) - [zeros(3);D_Y_X(R(t))]
% Compute <D_Y_X,Y>
g_DYX_Y = rotBundle_metric_nonNatural(Zt(t),D_Y_X(R(t)),Yt(t),m_nonnatural);

% Compute D_X_Y
D_Y_Y = rotBundle_covar_nonNatural(Y,Y,U_R,m_nonnatural);
D_Y_Y2 = rotRef_covar(@(R,U,RRef) [zeros(3);Y(R,U)],...
    @(R,U,RRef) [zeros(3);Y(R,U)],...
    R,U_R,RRef,M_test,t);
fprintf('Test rotRef_covar (should be zero), D_Y_Y at (R,U):\n');
D_Y_Y2(R(t)) - [zeros(3);D_Y_Y(R(t))]
% Compute <D_Y_Y,X>
g_DYY_X = rotBundle_metric_nonNatural(Zt(t),D_Y_Y(R(t)),Xt(t),m_nonnatural);

g_analytical_der = g_DYX_Y + g_DYY_X;
%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('analytical der: %0.5f:\n',g_analytical_der);
fprintf('Error of der: %0.5f:\n',g_nonnatural_der-g_analytical_der);