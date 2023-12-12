% Test if the change of the "non-natural" metric of two left-invarient 
% vector fields on  TSO(3)xSO(3) at (R,U,RRef) can be computed by using 
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
U = @(R,t) R*hat3(uVec);
du = 0*uVec;
% Define curve on reference manifold SO(3)
[RRef,~,RRef0,dRRef0] = rot_randGeodFun;
dRRef = @(RRef) RRef*RRef0'*dRRef0;
% Make m result in pos. def. matrix
M_nonnatural = randn(3,3);
M_nonnatural = M_nonnatural'*M_nonnatural; 
% Generate system dynamic vector field
x1 = randn(3,1);x2 = randn(3,1); x3 = randn(3,1);
X = @(R,U,RRef) [RRef*hat3(x1);R*hat3(x2);R*hat3(x3)];
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

% Compute D_Y_X in the "natural" coordinates
% IE D_Y_X_nonnatural = D_Y_X_natural + D_Y_Xm_natural + D_Ym_X_natural + D_Ym_Xm_natural
% where Xm and Ym are the additional VF on TSO(3) (ie Xm=J*X-X and Ym=J*Y-Y)
% Compute D_Y_X_natural at (R,RRef)
D_Y_X_natural = rotRef_covar(Y,X,R(t),U_R,RRef,M_natural,t);
% Compute D_Y_Xm_natural at (RRef,~)
Xm = @(RRef,U,R) [zeros(3); J(2,1)*RRef*hat3(x1); J(3,1)*RRef*hat3(x1)];
D_Y_Xm_natural = rotRef_covar(Y,Xm,RRef(t),U_R,RRef,M_natural,t);
% Compute D_Ym_X_natural at (R,~)
Ym = @(RRef,U,R) [zeros(3); J(2,1)*dRRef(RRef); J(3,1)*dRRef(RRef)];
D_Ym_X_natural = rotRef_covar(Ym,X,R(t),U_R,R,M_natural,t);
% Compute D_Ym_Xm_natural at (RRef,~)
D_Ym_Xm_natural = rotRef_covar(Ym,Xm,RRef(t),U_R,RRef,M_natural,t);
% Sum all partial covar (with translation to (R,RRef))
D_Y_X = @(R,RRef) D_Y_X_natural(R,RRef) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Y_Xm_natural(RRef,RRef) + ...
    D_Ym_X_natural(R,R) + LT_R(R,RRef)*LT_R(RRef,RRef)'*D_Ym_Xm_natural(RRef,RRef);

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
    
g_DYY_X_natural = rotRef_metric(Zt(t), ...
    D_Y_Y(R(t),RRef(t)), ...
    X(R(t),U(R(t),t),RRef(t)) + Xm(R(t),U(R(t),t),RRef(t)), ...
    M_natural);

g_natural_der = g_DYX_Y_natural + g_DYY_X_natural;

%% Compare the numerical approx der of d/dt<X,Y>_nonnatural to <D_Y_X,Y>
fprintf('numerical der: %0.5f:\n',g_nonnatural_der);
fprintf('analytical der: %0.5f:\n',g_natural_der);
fprintf('Error of der: %0.5f:\n',g_nonnatural_der-g_natural_der);