% Test if the change of the "non-natural" metric of a left-invarient 
% vector field and controlled system VF on  TSO(3)xSO(3) at (R,U,RRef) can be computed by using 
% left translation to the same tangent space for cross terms

% THIS DOES NOT WORK, see POC_rotRef_Sys_evalAtRRef.m

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
% Make m result in pos. def. matrix
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural; 
% Generate system dynamic vector field
X = @(R,U,RRef) [rot_log(RRef,eye(3));U;kd*rot_log(R,RRef)-kv*U];
% Generate vector field along curve Z(t)=[RRef;R;U];
Z = @(R,U,RRef) [RRef;R;U];
Y = @(R,U,RRef) [dRRef(RRef);dR(R);R*hat3(du)];
U_R = @(R) U(R,t); % Redefine U as only a function of R
LT_R = @(R,RRef) [RRef zeros(3,6);zeros(3) R zeros(3);zeros(3,6) R]; % Convient matrix for Left translation

% Find the natural metric and cooresponding transformation matrix J such
% that M_nonnatural = J'*M_natural*J
[J,M_natural] = rotRef_SchurComplement(M_nonnatural);
% M_nonnatural-J'*M_natural*J
%% Compute the change of the nonnatural metric (ie d/dt<X,Y>)
Xt = @(t) X(R(t),U(R(t),t),RRef(t));
Yt = @(t) Y(R(t),U(R(t),t),RRef(t));
Zt = @(t) Z(R(t),U(R(t),t),RRef(t));

g_nonnatural = @(t) rotRef_metric_nonNatural(Zt(t),Xt(t),Yt(t),M_nonnatural);
g_nonnatural_der = funApproxDer(g_nonnatural,t);

%% Compute the change of the nonnatural metric using covar ders (ie <D_Y_X,Y> + <D_Y_Y,X>)
% Define useful constants
R_t = R(t); U_t = U(R_t,t); RRef_t = RRef(t); % Evaluated at t
LT_TSO3 = @(R) kron(eye(2),R);

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

% Compute <D_Ym_Xm,Ym>
D_Ym_Xm_SO3 = rot_covar(Ym_SO3,Xm_SO3);
D_Ym_Xm_TSO3 = rotBundle_covar_nonNatural(Ym_TSO3_RU,Xm_TSO3_RU,U_R,m_natural_TSO3,'rigidrot',kv);
% Since Xm has dependencies on RRef, compute dXm/dRRef
dXm_dRRef = funApproxDer(@(t) Xm_TSO3(R_t,U_t,RRef(t)),t);
D_Ym_Xm_TSO3 = @(R) D_Ym_Xm_TSO3(R) + dXm_dRRef;
g_DYmXm_Ym_SO3 = M_natural(1,1)*rot_metric(RRef_t,D_Ym_Xm_SO3(RRef_t),Ym_SO3(RRef_t));
g_DYmXm_Ym_TSO3 = rotBundle_metric_nonNatural([R_t;U_t],...
    D_Ym_Xm_TSO3(R_t), Ym_TSO3(R_t,U_t,RRef_t),m_natural_TSO3);
g_DYX_Y_natural = g_DYmXm_Ym_SO3 + g_DYmXm_Ym_TSO3;

% Compute <D_Ym_Ym,Xm>
D_Ym_Ym_SO3 = rot_covar(Ym_SO3,Ym_SO3);
D_Ym_Ym_TSO3 = rotBundle_covar_nonNatural(Ym_TSO3_RU,Ym_TSO3_RU,U_R,m_natural_TSO3);
% Since Ym has dependencies on RRef, compute dYm/dRRef
dYm_dRRef = funApproxDer(@(t) Ym_TSO3(R_t,U_t,RRef(t)),t);
D_Ym_Ym_TSO3 = @(R) D_Ym_Ym_TSO3(R) + dYm_dRRef;
g_DYmYm_Xm_SO3 = M_natural(1,1)*rot_metric(RRef_t,D_Ym_Ym_SO3(RRef_t),Xm_SO3(RRef_t));
g_DYmYm_Xm_TSO3 = rotBundle_metric_nonNatural([R_t;U_t],...
    D_Ym_Ym_TSO3(R_t), Xm_TSO3(R_t,U_t,RRef_t),m_natural_TSO3);
g_DYY_X_natural = g_DYmYm_Xm_SO3 + g_DYmYm_Xm_TSO3;


% Sum the components to compute d/dt<X,Y>
g_natural = @(t) rotRef_metric(Zt(t),Xm(R(t),U(R(t),t),RRef(t)),...
    Ym(R(t),U(R(t),t),RRef(t)),M_natural);
g_natural_der = g_DYX_Y_natural + g_DYY_X_natural;

%% Compute the natural metric and its derivative separately
g_natural_SO3 = @(t) M_natural(1,1)*rot_metric(RRef(t),Xm_SO3(RRef(t)),Ym_SO3(RRef(t)));
g_natural_TSO3 = @(t) rotBundle_metric_nonNatural([R(t);U(R(t),t)],...
    Xm_TSO3(R(t),U(R(t),t),RRef(t)), Ym_TSO3(R(t),U(R(t),t),RRef(t)),...
    m_natural_TSO3);

g_natural_sep = @(t) g_natural_SO3(t) + g_natural_TSO3(t);
g_natural_SO3_der = funApproxDer(g_natural_SO3,t);
g_natural_TSO3_der = funApproxDer(g_natural_TSO3,t);
g_natural_sep_der = g_natural_SO3_der + g_natural_TSO3_der;

%% Compute the d/dt<X,Y> by decomposing <D_Y_X,Y> (only TSO3) into 2 terms by using linear property on 1st arg (similarly <D_Y_Y,X>)

% Splite Ym_TSO3 into two terms
Ym_TSO3_Y = @(R,U,RRef) [dR(R);R*hat3(du)];
Ym_TSO3_Y_RU = @(R,U) Ym_TSO3_Y(R,U,RRef(t));
Ym_TSO3_Yref = @(R,U,RRef) [J(2,1)*R*RRef'*dRRef(RRef);J(3,1)*R*RRef'*dRRef(RRef)];
Ym_TSO3_Yref_RU = @(R,U) Ym_TSO3_Yref(R,U,RRef(t));

% Compute <D_Y_Xm,Ym>
D_Y_Xm_TSO3 = rotBundle_covar_nonNatural(Ym_TSO3_Y_RU,Xm_TSO3_RU,U_R,m_natural_TSO3,'rigidrot',kv);
% Since Xm has dependencies on RRef, add additional term
D_Y_Xm_TSO3 = @(R)  D_Y_Xm_TSO3(R) + dXm_dRRef;
g_DYXm_Ym = rotBundle_metric_nonNatural([R(t);U(R(t),t)],...
    D_Y_Xm_TSO3(R(t)), Ym_TSO3(R(t),U(R(t),t),RRef(t)),...
    m_natural_TSO3);
% Compute <D_Yref_Xm,Ym>
% Assume we do not need to compute changes along fibers since Yref is
% really just changes along SO3 (no movement along fibers)
D_Yref_Xm_TSO3 = rotBundle_covar_nonNatural(Ym_TSO3_Yref_RU,Xm_TSO3_RU,U_R,m_natural_TSO3);
g_DYrefXm_Ym = rotBundle_metric_nonNatural([R(t);U(R(t),t)],...
    D_Yref_Xm_TSO3(R(t)), Ym_TSO3(R(t),U(R(t),t),RRef(t)),...
    m_natural_TSO3);

% Recompute <D_Y_X,Y> + <D_Y_Y,X>
g_DYX_Y_natural_2 = g_DYmXm_Ym_SO3 + g_DYXm_Ym + g_DYrefXm_Ym;
g_natural_der_2 = g_DYX_Y_natural_2 + g_DYY_X_natural;

%% Compute d/dt<X,Y> by spliting into covars evaluated at (R,U) and (RRef,0)

% Define "natural" vector fields on rotRef at (R,U,RRef) with left
% translation
Jm = J; Jm(1,1) = 0;Jm(2,2) = 0; Jm(3,3) = 0;
% Left translate vector on SO3 to R on TSO3
Xm_at_R_U_RRef = @(R,U,RRef) X(R,U,RRef) + LT_R(R,zeros(3))*LT_R(RRef,zeros(3))'*kron(Jm,eye(3))*X(R,U,RRef);
Ym_at_R_U_RRef = @(R,U,RRef) Y(R,U,RRef) + LT_R(R,zeros(3))*LT_R(RRef,zeros(3))'*kron(Jm,eye(3))*Y(R,U,RRef);
% Left translate vector on TSO3 to RRef on TSO3 (note: RRef1 = RRef2 is the
% original RRef)
Z_at_RRef_URRef_R = [RRef_t;RRef_t;zeros(3)];
Xm_at_RRef_URRef_R_funt = @(R,U,RRef) extractComp(LT_R(RRef,zeros(3))*LT_R(R,zeros(3))'*X(R,U,RRef)+...
    kron(Jm,eye(3))*X(R,U,RRef),4,9,1,3); % used to compute chain rule
Xm_at_RRef_URRef_R = @(RRef2,URRef2,RRef1) LT_R(RRef2,zeros(3))*LT_R(R(t),zeros(3))'*X(R(t),U(R(t),t),RRef2)+...
    kron(Jm,eye(3))*X(R(t),U(R(t),t),RRef2); % used to compute covar
Ym_at_RRef_URRef_R_funt = @(R,U,RRef) extractComp(LT_R(RRef,zeros(3))*LT_R(R,zeros(3))'*Y(R,U,RRef)+...
    kron(Jm,eye(3))*Y(R,U,RRef),4,9,1,3); % used to compute chain rule
Ym_at_RRef_URRef_R = @(RRef2,URRef2,RRef1) LT_R(RRef2,zeros(3))*LT_R(R(t),zeros(3))'*Y(R(t),U(R(t),t),RRef2)+...
    kron(Jm,eye(3))*Y(R(t),U(R(t),t),RRef2); % used to compute covar

% Define the Y vector as a sum that we'll use to splite the covar into two
% parts in the first argument (IE D_Y_X = D_Y1_X|(R,U) + D_Y2_X|(RRef,0)
% where Y1 = Y at (R,U,RRef) and Y2 = kron(Jm,eye(3))*Y(R,U,RRef) at
% (RRef,URRef,R) (NOTE: Y1+Y2 (with left translation) is Ym)
Y1 = Y; % Vector in the TSO3 slot is at (R,U)
Y2 = @(RRef2,URRef2,RRef1) kron(Jm,eye(3))*Y(R(t),U(R(t),t),RRef2); % Vector in the TSO3 slot is at (RReef,URRef)

% Compute <D_Y_X,Y>_natural = <D_Y1_X1,Ym_at_R_U_RRef> + <D_Y2_X2,Ym_at_R_U_RRef>
% where X1 = Xm_at_R_U_RRef and X2 = Xm_at_RRef_URRef_R
% <D_Y1_X1,Ym_at_R_U_RRef>
D_Y1_X1 = rotRef_covar(Y1,Xm_at_R_U_RRef,R,U_R,RRef,M_natural,t,'rigidrot',kv);
g_DY1X1_Ym = rotRef_metric(Zt(t),D_Y1_X1(R_t),Ym_at_R_U_RRef(R_t,U_t,RRef_t),M_natural);
% <D_Y2_X2,Ym_at_R_U_RRef>
D_Y2_X2 = rotRef_covar(Y2,Xm_at_RRef_URRef_R,R,U_R,RRef,M_natural,t,...
    'evalatrrefontso3',@(t) Xm_at_RRef_URRef_R_funt(R(t),U(R(t),t),RRef_t));
g_DY2X2_Ym = rotRef_metric(Z_at_RRef_URRef_R,D_Y2_X2(RRef_t),...
    Ym_at_RRef_URRef_R(RRef_t,zeros(3),RRef_t),M_natural);
% Compute <D_Y_X,Y>_natural
g_DYX_Y_new = g_DY1X1_Ym + g_DY2X2_Ym;

% Compute <D_Y_Y,X>_natural = <D_Y1_Y1,Xm_at_R_U_RRef> + <D_Y2_Y2,Xm_at_RRef_URRef_R> 
% <D_Y1_Ym,Xm_at_R_U_RRef>
D_Y1_Y1 = rotRef_covar(Y1,Y1,R,U_R,RRef,M_natural,t);
g_DY1Y1_Xm = rotRef_metric(Zt(t),D_Y1_Y1(R_t),Xm_at_R_U_RRef(R_t,U_t,RRef_t),M_natural);
% <D_Y2_Y2,Xm_at_RRef_URRef_R> 
D_Y2_Y2 = rotRef_covar(Y2,Y2,R,U_R,RRef,M_natural,t,...
    'evalatrrefontso3', @(t) Ym_at_RRef_URRef_R_funt(R(t),U(R(t),t),RRef_t));
g_DY2Y2_Xm = rotRef_metric(Z_at_RRef_URRef_R,D_Y2_Y2(RRef_t),...
    Xm_at_RRef_URRef_R(RRef_t,U_t,RRef_t),M_natural);
% Compute <D_Y_Y,X>_natural
g_DYY_X_new = g_DY1Y1_Xm + g_DY2Y2_Xm;

% Compute d/dt<X,Y>_natural
g_natural_der_new = g_DYX_Y_new+g_DYY_X_new

%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('Nonnatural metric: %0.5f:\n',g_nonnatural(t));
fprintf('Natural metric: %0.5f:\n',g_natural(t));
fprintf('Natural metric (sep): %0.5f:\n',g_natural_sep(t));
fprintf('Error of metrics: %0.5f:\n',g_nonnatural(t)-g_natural(t));
fprintf('Numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('Analytical der: %0.5f:\n',g_natural_der);
fprintf('Analytical der (sep): %0.5f:\n',g_natural_der_2);
fprintf('Numerical der (sep): %0.5f:\n',g_natural_sep_der);
fprintf('Error of ders: %0.5f:\n',g_nonnatural_der-g_natural_der);

%% Compare the natural metric derivative analytically and with numerical results
fprintf('Error der SO3: %0.5f:\n', g_natural_SO3_der - g_DYmXm_Ym_SO3+g_DYmYm_Xm_SO3);
fprintf('Error der TSO3: %0.5f:\n', g_natural_TSO3_der - g_DYmYm_Xm_TSO3-g_DYmXm_Ym_TSO3);

