% Test vector fields on TSO3 with a componenet dependent on a reference SO3
% manifold (by left translation), consider the case where the first
% argument is only the transformed vector from Ref

close all; clear all; clc;
fprintf('Test Script: "%s"\n\n',mfilename)
%% Define Parameters
t = rand; % Pick a random time to evaluate
% Define a stationary curve on TSO(3)
% R0 = rot_randn;
% R = @(t) R0;
% dR = @(R) zeros(3);
% % Linear U
% uVec = randn(3,1);
% U = @(R,t) R*hat3(uVec);
% du = 0*uVec;
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
X = @(R,U,RRef) [U+R*RRef'*rot_log(RRef,eye(3));5*rot_log(R,RRef)-2*U+R*RRef'*rot_log(RRef,eye(3))];
X_LT_RRef = @(R,U,RRef) [RRef*R'*U+rot_log(RRef,eye(3));5*RRef*R'*rot_log(R,RRef)-2*RRef*R'*U+rot_log(RRef,eye(3))]; % Left translate X to (RRef)
% X = @(R,U,RRef) [U;R*hat3(x1)+R*hat3(x2)];
X_R = @(R,U) X(R,U,RRef(t));
X_RRef = @(RRef,URef) [RRef*R(t)'*U(R(t),t)+rot_log(RRef,eye(3));5*RRef*R(t)'*rot_log(R(t),RRef)-2*RRef*R(t)'*U(R(t),t)+rot_log(RRef,eye(3))]; % for use when computing covar at RRef on TSO3
% Generate vector field along the curve Z(t) = [R(t);U(R,t)]
Z = @(R,U) [R;U];
Y = @(R,U,RRef) [dR(R)+R*RRef'*dRRef(RRef);R*hat3(du)+R*RRef'*dRRef(RRef)]; % NOTE: dZ/dt = Y;
Y_R = @(R,U) Y(R,U,RRef(t));
Y_RRef = @(RRef,URef) [dRRef(RRef);dRRef(RRef)]; % for use when computing covar at RRef on TSO3
Y_RRef_all = @(R,U,RRef) [RRef*R'*dR(R)+dRRef(RRef);RRef*R'*R*hat3(du)+dRRef(RRef)];
Y_RRef_all_RRef = @(RRef,URef) Y_RRef_all(R(t),U(R(t),t),RRef);
U_R = @(R) U(R,t);
LT_R = @(R) kron(eye(2),R); % Convient matrix for Left translation

%% Compute the change of the nonnatural metric (ie d/dt<X,Y>)
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t),RRef(t));
Zt = @(t) Z(R(t),U(R(t),t));

g_nonnatural = @(t) rotBundle_metric_nonNatural(Zt(t),Xt(t),Yt(t),m_nonnatural);
g_nonnatural_der = funApproxDer(g_nonnatural,t);
%% Compute the change of the nonnatural metric using covar ders (ie <D_Y_X,Y> + <D_Y_Y,X>)
% Define useful constants
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t); % Evaluated at t

% Compute D_Y_X (note X_RRef does not depend on U_RRef so no need to
% account for changes along fibers at RRef)
D_Y_X = rotBundle_covar_nonNatural(Y_RRef,X_RRef,@(RRef) zeros(3),m_nonnatural); % assume URef = 0 since X is consider "constant" when only moving on Ref
% Additional term to account for the tangent vector from R depending on
% R (this is seen as a direct time dependence on TSO3)
dX_dt = funApproxDer(@(t) X_LT_RRef(R(t),U(R(t),t),RRef_t),t);
D_Y_X = @(R) D_Y_X(R) + dX_dt;
M_test = zeros(3,3); M_test(2:3,2:3)=M_nonnatural;
D_Y_X2 = rotRef_covar(@(RRef2,URRef2,RRef1) [zeros(3);Y_RRef(RRef2,URRef2)],...
    @(RRef2,URRef2,RRef1) [zeros(3);X_RRef(RRef2,URRef2)],...
    R, U, RRef, M_test,t,'evalatrrefontso3',@(t) X_LT_RRef(R(t),U(R(t),t),RRef_t));
fprintf('Test rotRef_covar (should be zero), D_Y_X at (RRef,0):\n');
D_Y_X2(RRef(t)) - [zeros(3);D_Y_X(RRef(t))]
% Compute <D_Y_X,Y> at (RRef,0)
% g_DYX_Y = rotBundle_metric_nonNatural([RRef(t);zeros(3)],D_Y_X(RRef(t)),LT_R(RRef(t))*LT_R(R(t))'*Yt(t),m_nonnatural);
g_DYX_Y = rotBundle_metric_nonNatural([RRef(t);zeros(3)],D_Y_X(RRef(t)),Y_RRef_all_RRef(RRef_t,zeros(3)),m_nonnatural);
g_DYX_Y2 = rotRef_metric([RRef(t);RRef(t);zeros(3)],D_Y_X2(RRef(t)),[zeros(3);LT_R(RRef(t))*LT_R(R(t))'*Yt(t)],M_test);

% Compute D_X_Y
D_Y_Y = rotBundle_covar_nonNatural(Y_RRef,Y_RRef_all_RRef,@(RRef) zeros(3),m_nonnatural);
dYall_dt = funApproxDer(@(t) Y_RRef_all(R(t),U(R(t),t),RRef_t),t);
D_Y_Y = @(R) D_Y_Y(R) + dYall_dt;
D_Y_Y2 = rotRef_covar(@(RRef2,URRef2,RRef1) [zeros(3);Y_RRef(RRef2,URRef2)],...
    @(RRef2,URRef2,RRef1) [zeros(3);Y_RRef_all_RRef(RRef2,URRef2)],...
    R, U, RRef, M_test,t,'evalatrrefontso3',@(t) Y_RRef_all(R(t),U(R(t),t),RRef_t)   );
fprintf('Test rotRef_covar (should be zero), D_Y_Y at (RRef,0):\n');
D_Y_Y2(RRef(t)) - [zeros(3);D_Y_Y(RRef(t))]
% Compute <D_Y_Y,X> at (RRef,0)
% g_DYY_X = rotBundle_metric_nonNatural([RRef(t);zeros(3)],D_Y_Y(RRef(t)),LT_R(RRef(t))*LT_R(R(t))'*Xt(t),m_nonnatural);
g_DYY_X = rotBundle_metric_nonNatural([RRef(t);zeros(3)],D_Y_Y(RRef(t)),X_LT_RRef(R_t,U_t,RRef_t),m_nonnatural);

g_analytical_der = g_DYX_Y + g_DYY_X;
%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('analytical der: %0.5f:\n',g_analytical_der);
fprintf('Error of der: %0.5f:\n',g_nonnatural_der-g_analytical_der);