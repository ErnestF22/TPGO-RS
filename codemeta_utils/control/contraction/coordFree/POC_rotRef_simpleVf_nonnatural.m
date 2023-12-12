% Test if the change of the "non-natural" metric of a left-invarient 
% vector field and "simple" VF on  TSO(3)xSO(3) at (R,U,RRef) can be computed by using 
% left translation to the same tangent space for cross terms

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
% Make m result in pos. def. matrix
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural; 
% Generate system dynamic vector field
x1 = randn(3,1);x2 = randn(3,1);
X = @(R,U,RRef) [RRef*hat3(x1);U;R*hat3(x2)-2*U];
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
R0 = R(t); UR0 = U(R0,t); RRef0 = RRef(t); URRef0 = U(RRef0,t);

% Compute D_Y_X in the "natural" coordinates
% IE D_Y_X_nonnatural = D_Y_X_natural + D_Y_Xm_natural + D_Ym_X_natural + D_Ym_Xm_natural
% where Xm and Ym are the additional VF on TSO(3) (ie Xm=J*X-X and Ym=J*Y-Y)
% Compute D_Y_X_natural at (R,RRef)
D_Y_X_natural = rotRef_covar(Y,X,R(t),U_R,RRef,M_natural,t,'rigidrot',2);
% Compute D_Y_Xm_natural at (RRef,~)
Xm = @(R,U,RRef) [zeros(3); J(2,1)*R*hat3(x1); J(3,1)*R*hat3(x1)];
D_Y_Xm_natural = rotRef_covar(Y,Xm,RRef(t),U_R,RRef,M_natural,t);
% Compute D_Ym_X_natural at (R,~)
Ym = @(R,U,RRef) [zeros(3); J(2,1)*dRRef(R); J(3,1)*dRRef(R)];
D_Ym_X_natural = rotRef_covar(Ym,X,R(t),U_R,R,M_natural,t);
% Compute D_Ym_Xm_natural at (RRef,~)
D_Ym_Xm_natural = rotRef_covar(Ym,Xm,RRef(t),U_R,RRef,M_natural,t);
% Sum all partial covar (with translation to (R,RRef))
D_Y_X = @(R,RRef) D_Y_X_natural(R,RRef) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Y_Xm_natural(RRef,RRef) + ...
    D_Ym_X_natural(R,R) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Ym_Xm_natural(RRef,RRef);

% Test linearlity implemention of rotRef covar
Ynew = @(R,U,RRef) Y(R,U,RRef)+Ym(R,U,RRef);
D_Ynew_X = rotRef_covar(Ynew,X,R(t),U_R,RRef,M_natural,t,'rigidrot',2);
SumCovar_YplusYm_X = @(R,RRef) D_Y_X_natural(R,RRef) + D_Ym_X_natural(R,RRef);
fprintf('(linear in 1st arg) Should be zero:\n')
D_Ynew_X(R(t),RRef(t))-SumCovar_YplusYm_X(R(t),RRef(t))
%Above works, but not "correct"
% Do we really ignore the covar along the fibers in the direction of
% J(3,1)*dRRef(R) since J(3,1)*dRRef(R) is only movement about the
% reference SO3 manifold?

% Test addititive prop of covar of 2nd arg
Xnew = @(R,U,RRef) X(R,U,RRef) + Xm(R,U,RRef);
D_Y_Xnew = rotRef_covar(Y,Xnew,R(t),U_R,RRef,M_natural,t,'rigidrot',2);
SumCovar_Y_XplusXm = @(R,RRef) D_Y_X_natural(R,RRef) + D_Y_Xm_natural(R,RRef);
fprintf('(additive in 2nd arg) Should be zero:\n')
D_Y_Xnew(R(t),RRef(t)) - SumCovar_Y_XplusXm(R(t),RRef(t))
% Above works

% Compute D_Y_Y in the "natural" coordinates
% Compute D_Y_Y_natural at (R,RRef)
D_Y_Y_natural = rotRef_covar(Y,Y,R(t),U_R,RRef,M_natural,t);
% Compute D_Y_Ym_natural at (RRef,~)
D_Y_Ym_natural = rotRef_covar(Y,Ym,RRef(t),U_R,RRef,M_natural,t);
% Compute D_Ym_Y_natural at (R,~)
D_Ym_Y_natural = rotRef_covar(Ym,Y,R(t),U_R,R,M_natural,t);
% Compute D_Ym_Ym_natural at (RRef,~)
D_Ym_Ym_natural = rotRef_covar(Ym,Ym,RRef(t),U_R,RRef,M_natural,t);
% Sum all partial covar (with translation to (R,RRef))
D_Y_Y = @(R,RRef) D_Y_Y_natural(R,RRef) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Y_Ym_natural(RRef,RRef) + ...
    D_Ym_Y_natural(R,R) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Ym_Ym_natural(RRef,RRef);

% Compute <D_Y_X,Y> + <D_Y_Y,X> in "natural" coords
g_DYX_Y_natural = rotRef_metric(Zt(t), ...
    D_Y_X(R(t),RRef(t)),...
    Y(R(t),U(R(t),t),RRef(t)) + Ym(R(t),U(R(t),t),RRef(t)), ...
    M_natural);    
% Compute metric at their respective location (without translating)
g_DYX_natural = rotRef_metric(Zt(t),...
    D_Y_X_natural(R(t),RRef(t)),...
    Y(R(t),U(R(t),t),RRef(t)) + Ym(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
g_DYXm_natural = rotRef_metric([RRef(t);RRef(t);U(RRef(t),t)],...
    D_Y_Xm_natural(RRef(t),RRef(t)),...
    Y(RRef(t),U(RRef(t),t),RRef(t)) + Ym(RRef(t),U(RRef(t),t),RRef(t)), ...
    M_natural);
g_DYmX_natural = rotRef_metric([R(t);R(t);U(R(t),t)],...
    D_Ym_X_natural(R(t),R(t)),...
    Y(R(t),U(R(t),t),R(t)) + Ym(R(t),U(R(t),t),R(t)), ...
    M_natural);
g_DYmXm_natural = rotRef_metric([RRef(t);RRef(t);U(RRef(t),t)],...
    D_Ym_Xm_natural(RRef(t),RRef(t)),...
    Y(RRef(t),U(RRef(t),t),RRef(t)) + Ym(RRef(t),U(RRef(t),t),RRef(t)), ...
    M_natural);
g_DYX_Y_natural_diff_space = g_DYX_natural + g_DYXm_natural + g_DYmX_natural + g_DYmXm_natural;

% <D_Y_Y,X>
g_DYY_X_natural = rotRef_metric(Zt(t), ...
    D_Y_Y(R(t),RRef(t)), ...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
% Compute metric at (R,RRef) since X is not left-invar and D_Y?_Y? is
g_DYY_natural = rotRef_metric(Zt(t),...
    D_Y_Y_natural(R(t),RRef(t)),...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
g_DYYm_natural = rotRef_metric(Zt(t),...
    D_Y_Ym_natural(R(t),RRef(t)),...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
g_DYmY_natural = rotRef_metric(Zt(t),...
    D_Ym_Y_natural(R(t),RRef(t)),...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
g_DYmYm_natural = rotRef_metric(Zt(t),...
    D_Ym_Ym_natural(R(t),RRef(t)),...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);
g_DYY_X_natural_diff_space = g_DYY_natural + g_DYYm_natural + g_DYmY_natural + g_DYmYm_natural;

g_natural_der = g_DYX_Y_natural + g_DYY_X_natural;
g_nautral_der_diff_space = g_DYX_Y_natural_diff_space + g_DYY_X_natural_diff_space;
%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('analytical der: %0.5f:\n',g_natural_der);
fprintf('diff_space der: %0.5f:\n',g_nautral_der_diff_space);
fprintf('Error of der: %0.5f:\n',g_nonnatural_der-g_natural_der);