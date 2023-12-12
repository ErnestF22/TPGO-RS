% Test if the change of the "non-natural" metric of a left-invarient 
% vector field and controlled system VF on  TSO(3)xSO(3) at (R,U,RRef) can be computed by using 
% left translation to the same tangent space for cross terms

close all; clear all; clc;
fprintf('Test Script: "%s"\n\n',mfilename)
%% Define Parameters
t = rand; % Pick a random time to evaluate
kd = 5; kv = 2;
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
% RRef = @(t) RRef0; dRRef = @(R) zeros(3); % Uncomment to make curve on SO3 stationary
% Make m result in pos. def. matrix
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural; 
% Generate system dynamic vector field
X = @(R,U,RRef) [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];
% Generate vector field along curve Z(t)=[RRef;R;U];
Z = @(R,U,RRef) [RRef;R;U];
Y = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(du)];
% Y = @(R,U,RRef) [dRRef(RRef);zeros(6,3)]; R = @(t) R0; dR = @(R) zeros(3); U = @(R,t) R*hat3(uVec); du = 0*uVec; % Uncommeent to test only motion on SO3
U_R = @(R) U(R,t); % Redefine U as only a function of R
LT_R = @(R,RRef) [RRef zeros(3,6);zeros(3) R zeros(3);zeros(3,6) R]; % Convient matrix for Left translation

% Find the natural metric and cooresponding transformation matrix J such
% that M_nonnatural = J'*M_natural*J
[J,M_natural] = rotRef_SchurComplement(M_nonnatural);
% M_nonnatural-J'*M_natural*J
%% Compute the change of the nonnatural metric (ie d/dt<X,Y>)
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t); % Evaluated at t
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t),RRef(t));
Zt = @(t) Z(R(t),U(R(t),t),RRef(t));

g_nonnatural = @(t) rotRef_metric_nonNatural(Zt(t),Xt(t),Yt(t),M_nonnatural);
g_nonnatural_der = funApproxDer(g_nonnatural,t);

%% Compute the change of the nonnatural metric using covar ders (ie <D_Y_X,Y> + <D_Y_Y,X>)
% Define useful constants
LT_TSO3 = @(R) kron(eye(2),R);

% Compute d/dt<X,Y> for SO3 component
Xm_SO3 = @(RRef) rot_log(RRef,eye(3));
Ym_SO3 = @(RRef) dRRef(RRef);
D_Ym_Xm_SO3 = rot_covar(Ym_SO3,Xm_SO3);
g_DYmXm_Ym_SO3 = M_natural(1,1)*rot_metric(RRef_t,D_Ym_Xm_SO3(RRef_t),Ym_SO3(RRef_t));
D_Ym_Ym_SO3 = rot_covar(Ym_SO3,Ym_SO3);
g_DYmYm_Xm_SO3 = M_natural(1,1)*rot_metric(RRef_t,D_Ym_Ym_SO3(RRef_t),Xm_SO3(RRef_t));
g_natural_der_SO3 = g_DYmXm_Ym_SO3 + g_DYmYm_Xm_SO3;

%% Compute d/dt<X,Y> on TSO3 by splitting into evaluation at (R,U) and (RRef,URRef)
% IE d/dt<X,Y>_nonnatural = <D_Y_JX,JY>_(R,U) + <D_Y_JY,JX>_(R,U) +
% <D_Y'_JX,JY>_(RRef,0) + <D_Y'_JY,JX>_(RRef,0) where JY = Y+Y'

% Compute at (R,U)
Zt_TSO3_R = [R(t);U(R(t),t)];
m_nonnatural = [M_natural(2,2);M_natural(2,3);M_natural(3,3)];
X_TSO3_R_all = @(R,U,RRef) [U+J(2,1)*R*RRef'*rot_log(RRef,eye(3));kd*rot_log(R,RRef)-kv*U+J(3,1)*R*RRef'*rot_log(RRef,eye(3))];
X_TSO3_R_arg2 = @(R,U) X_TSO3_R_all(R,U,RRef_t);
Y_TSO3_R_all = @(R,U,RRef) [dR(R)+J(2,1)*R*RRef'*dRRef(RRef);R*hat3(du)+J(3,1)*R*RRef'*dRRef(RRef)];
Y_TSO3_R_arg1 = @(R,U) Y_TSO3_R_all(R,U,zeros(3));
Y_TSO3_R_arg2 = @(R,U) Y_TSO3_R_all(R,U,RRef_t);
% <D_Y_JX,JY>
D_Y_JX_R = rotBundle_covar_nonNatural(Y_TSO3_R_arg1,X_TSO3_R_arg2,U_R,m_nonnatural,'rigidrot',kv);
% Chain rule
dJX_dt_R = funApproxDer(@(t) X_TSO3_R_all(R_t,U_t,RRef(t)),t);
D_Y_JX_R = @(R) D_Y_JX_R(R) + dJX_dt_R;
g_DYJX_JY = rotBundle_metric_nonNatural(Zt_TSO3_R,D_Y_JX_R(R_t),Y_TSO3_R_all(R_t,U_t,RRef_t),m_nonnatural);
% <D_Y_JY,JX>
D_Y_JY_R = rotBundle_covar_nonNatural(Y_TSO3_R_arg1,Y_TSO3_R_arg2,U_R,m_nonnatural);
% Chain rule
dJY_dt_R = funApproxDer(@(t) Y_TSO3_R_all(R_t,U_t,RRef(t)),t);
D_Y_JY_R = @(R) D_Y_JY_R(R) + dJY_dt_R;
g_DYJY_JX = rotBundle_metric_nonNatural(Zt_TSO3_R,D_Y_JY_R(R_t),X_TSO3_R_all(R_t,U_t,RRef_t),m_nonnatural);
% Sum terms to get covar at (R,U) when 1st arg is only the component
% intially on TSO3
g_natural_TSO3_R = g_DYJX_JY + g_DYJY_JX;

% Try to reconstruct above results using functions
% Extract the individual vectors
X_RRef_fun = @(R,U,RRef) extractComp(X(R,U,RRef),1,3,1,3); % \in T_{RRef}SO3
X_R_fun = @(R,U,RRef) extractComp(X(R,U,RRef),4,6,1,3); % \in T_{R}SO3
X_U_fun = @(R,U,RRef) extractComp(X(R,U,RRef),7,9,1,3); % \in T_{R}SO3
Y_RRef_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),1,3,1,3); % \in T_{RRef}SO3
Y_R_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),4,6,1,3); % \in T_{R}SO3
Y_U_fun = @(R,U,RRef) extractComp(Y(R,U,RRef),7,9,1,3); % \in T_{R}SO3
% Construct the coordinate changed vectors
X_TSO3_R_all_fun = @(R,U,RRef) [X_R_fun(R,U,RRef) + J(2,1)*R*RRef'*X_RRef_fun(R,U,RRef);...
    X_U_fun(R,U,RRef) + J(3,1)*R*RRef'*X_RRef_fun(R,U,RRef)];
X_TSO3_R_arg2_fun = @(R,U) X_TSO3_R_all_fun(R,U,RRef_t);
Y_TSO3_R_all_fun = @(R,U,RRef) [Y_R_fun(R,U,RRef) + J(2,1)*R*RRef'*Y_RRef_fun(R,U,RRef);...
    Y_U_fun(R,U,RRef) + J(3,1)*R*RRef'*Y_RRef_fun(R,U,RRef)];
Y_TSO3_R_arg1_fun = @(R,U) Y_TSO3_R_all_fun(R,U,zeros(3));
Y_TSO3_R_arg2_fun = @(R,U) Y_TSO3_R_all_fun(R,U,RRef_t);
% <D_Y_JX,JY>
D_Y_JX_R_fun = rotBundle_covar_nonNatural(Y_TSO3_R_arg1_fun, X_TSO3_R_arg2_fun,U_R,m_nonnatural,'rigidrot',kv);
% Chain rule wrt RRef(t)
dJX_dt_R_fun = funApproxDer(@(t) X_TSO3_R_all_fun(R_t,U_t,RRef(t)),t);
D_Y_JX_R_fun = @(R) D_Y_JX_R_fun(R) + dJX_dt_R_fun;
g_DYJX_JY_fun = rotBundle_metric_nonNatural(Zt_TSO3_R,D_Y_JX_R_fun(R_t),Y_TSO3_R_all_fun(R_t,U_t,RRef_t),m_nonnatural);
% <D_Y_JY,JX>
D_Y_JY_R_fun = rotBundle_covar_nonNatural(Y_TSO3_R_arg1_fun,Y_TSO3_R_arg2_fun,U_R,m_nonnatural);
% Chain rule wrt RRef(t)
dJY_dt_R_fun = funApproxDer(@(t) Y_TSO3_R_all_fun(R_t,U_t,RRef(t)),t);
D_Y_JY_R_fun = @(R) D_Y_JY_R_fun(R) + dJY_dt_R_fun;
g_DYJY_JX_fun = rotBundle_metric_nonNatural(Zt_TSO3_R,D_Y_JY_R_fun(R_t),X_TSO3_R_all_fun(R_t,U_t,RRef_t),m_nonnatural);
% Sum terms to get covar at (R,U) when 1st arg is only the component
% intially on TSO3
g_natural_TSO3_R_fun = g_DYJX_JY_fun + g_DYJY_JX_fun;
% Compute the total der using the functions
g_natural_total_der_fun = g_natural_der_SO3 + g_natural_TSO3_R_fun;
fprintf('Analytical der(using functions): %0.5f\n\n', g_natural_total_der_fun);

%%%%%%%RESULTS INDICATE BELOW SHOULD BE 0, REASON: the curve doesn't change
%%%%%%%under coordinate transformation, thus the first argument of the
%%%%%%%covar metric compatibility test is unchanged. However, the Y vector
%%%%%%%field as the 2nd covar arg and metric arg are both transformed by
%%%%%%%the coordinate transformation.
% Compute at (RRef,0)
Zt_TSO3_RRef = [RRef_t;zeros(3)];
X_TSO3_RRef_all = @(R,U,RRef) [RRef*R'*U+J(2,1)*rot_log(RRef,eye(3));kd*RRef*R'*rot_log(R,RRef)-kv*RRef*R'*U+J(3,1)*rot_log(RRef,eye(3))];
X_TSO3_RRef_arg2 = @(RRef,URef) X_TSO3_R_all(R_t,U_t,RRef);
Y_TSO3_RRef_all = @(R,U,RRef) [RRef*R'*dR(R)+J(2,1)*dRRef(RRef);RRef*R'*R*hat3(du)+J(3,1)*dRRef(RRef)];
Y_TSO3_RRef_arg1 = @(RRef,URef) Y_TSO3_RRef_all(zeros(3),zeros(3),RRef);
Y_TSO3_RRef_arg2 = @(RRef,URef) Y_TSO3_RRef_all(R_t,U_t,RRef);
% <D_Y'_JX,JY>
D_Yprime_JX_RRef = rotBundle_covar_nonNatural(Y_TSO3_RRef_arg1,X_TSO3_RRef_arg2,@(RRef) zeros(3),m_nonnatural);
% Chain rule
dJX_dt_RRef = funApproxDer(@(t) X_TSO3_RRef_all(R(t),U(R(t),t),RRef_t),t);
D_Yprime_JX_RRef = @(R) D_Yprime_JX_RRef(R) + dJX_dt_RRef;
g_DYprimeJX_JY_RRef = rotBundle_metric_nonNatural(Zt_TSO3_RRef,D_Yprime_JX_RRef(RRef_t),Y_TSO3_RRef_all(R_t,U_t,RRef_t),m_nonnatural);
% <D_Y'_JY,JX>
D_Yprime_JY_RRef = rotBundle_covar_nonNatural(Y_TSO3_RRef_arg1,Y_TSO3_RRef_arg2,@(RRef) zeros(3),m_nonnatural);
% Chain rule
dJY_dt_RRef = funApproxDer(@(t) Y_TSO3_RRef_all(R(t),U(R(t),t),RRef_t),t);
D_Yprime_JY_RRef = @(R) D_Yprime_JY_RRef(R) + dJY_dt_RRef;
g_DYprimeJY_JX = rotBundle_metric_nonNatural(Zt_TSO3_RRef,D_Yprime_JY_RRef(RRef_t),X_TSO3_RRef_all(R_t,U_t,RRef_t),m_nonnatural);

fprintf('Covar with Yprime as first arg result is: %0.5f\n\n', g_DYprimeJX_JY_RRef + g_DYprimeJY_JX);
%%%%%%%RESULTS INDICATE ABOVE SHOULD BE 0
g_natural_TSO3_RRef = 0;

% Compute the total d/dt<X,Y>_natural
g_natural_der_TSO3 = g_natural_TSO3_R + g_natural_TSO3_RRef;
g_natural_total_der = g_natural_der_SO3 + g_natural_der_TSO3;

%% Compare results
fprintf('Numerical der: %0.5f\n', g_nonnatural_der);
fprintf('Analytical der: %0.5f\n', g_natural_total_der);
fprintf('Error der: %0.5f\n', g_nonnatural_der - g_natural_total_der);

%% Compute the natural metric and its derivative separately

% Define the transformed vector fields on their respective manifold
Xm_SO3 = @(RRef) rot_log(RRef,eye(3));
Xm_TSO3 = @(R,U,RRef) [U+J(2,1)*R*RRef'*rot_log(RRef,eye(3));...
    kd*rot_log(R,RRef)-kv*U+J(3,1)*R*RRef'*rot_log(RRef,eye(3))];
Xm_TSO3_RU = @(R,U) Xm_TSO3(R,U,RRef(t)); % Redefine assuming RRef = const
Xm = @(R,U,RRef) [Xm_SO3(RRef);Xm_TSO3(R,U,RRef)]; % 
Ym_SO3 = @(RRef) dRRef(RRef);
Ym_TSO3 = @(R,U,RRef) [dR(R)+J(2,1)*R*RRef'*dRRef(RRef);...
    R*hat3(du)+J(3,1)*R*RRef'*dRRef(RRef)];
Ym_TSO3_RU = @(R,U) Ym_TSO3(R,U,RRef(t)); % Redefine assuming RRef = const
Ym = @(R,U,RRef) [Ym_SO3(RRef);Ym_TSO3(R,U,RRef)];
m_natural_TSO3 = [M_natural(2,2);M_natural(2,3);M_natural(3,3)];
g_natural_SO3 = @(t) M_natural(1,1)*rot_metric(RRef(t),Xm_SO3(RRef(t)),Ym_SO3(RRef(t)));
g_natural_TSO3 = @(t) rotBundle_metric_nonNatural([R(t);U(R(t),t)],...
    Xm_TSO3(R(t),U(R(t),t),RRef(t)), Ym_TSO3(R(t),U(R(t),t),RRef(t)),...
    m_natural_TSO3);

g_natural_sep = @(t) g_natural_SO3(t) + g_natural_TSO3(t);
g_natural_SO3_der = funApproxDer(g_natural_SO3,t);
g_natural_TSO3_der = funApproxDer(g_natural_TSO3,t);
g_natural_sep_der = g_natural_SO3_der + g_natural_TSO3_der;

% Print results
fprintf('Numerical der (natural): %0.5f\n', g_natural_sep_der);
fprintf('Error der (SO3): %0.5f\n', g_natural_SO3_der - g_natural_der_SO3);
fprintf('Error der (TSO3): %0.5f\n', g_natural_TSO3_der - g_natural_der_TSO3);