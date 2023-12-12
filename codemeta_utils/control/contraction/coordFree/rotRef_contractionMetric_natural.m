function [g_DJYJX_JY] = rotRef_contractionMetric_natural(R,U,RRef,Y,M_nn,kd,kv)
% Compute the contraction metric using natural coordinates ( ie. <D_Y_X,Y>)
% INPUTS:
%   R, U, RRef := point on TSO(3)xSO(3) to evaluate
%   t := scalar time of evaluation
%   Y := The random vector to test [R*hat3(nu);R*hat3(zeta);R*hat3(eta)]
%       given in nonnatural coordinates
%   M_nn := the nonnatural metric [3x3] symmetric matrix
%   kd, kv := scalar positive control gains
% OUTPUTS:
%   g_DJYJX_JY := the metric computed at (R,U,RRef) in natural
%       coordinates

% Extract useful parameters
nu = rot_vee(RRef, extractComp(Y,1,3,1,3));
zeta = rot_vee(R, extractComp(Y,4,6,1,3));
eta = rot_vee(R, extractComp(Y,7,9,1,3));
w = rot_vee(R,U);
% Find the necessary coordinate transformation matrix
[J,M_n] = rotRef_SchurComplement(M_nn);
alpha = J(2,1); beta = J(3,1); gamma = M_n(1,1);
m1 = M_n(2,2); m2 = M_n(2,3); m3 = M_n(3,3);

% Compute gamma*<D_Y_X,Y>|_{SO3}
g_D_Y_X_SO3 = -gamma*nu'*rot3_logDiffMat(eye(3),RRef)*nu;
% Compute <m1*JYh+m2*JYv,D_JYh_JXh>
g_DJYh_JXh = -1/2*m2*(eta+beta*nu)'*(hat3(w)-alpha*rot_log(eye(3),RRef))*(zeta+alpha*nu)...
    -( m1*(zeta+alpha*nu)+m2*(eta+beta*nu) )'*alpha*rot3_logDiffMat(eye(3),RRef)*nu...
    +( m1*(zeta+alpha*nu)+m2*(eta+beta*nu) )'*(eta+beta*nu);
% Compute m2*<R(U,JYh)JXh,JYh>
g_R_U_JYh_JXh = m2/4*(zeta+alpha*nu)'*hat3(w)*(hat3(w)-alpha*rot_log(eye(3),RRef))*(zeta+alpha*nu);
% Compute -m3*<JYv,R(JYh,JXh)U>
g_R_JYh_JXh_U = m3/4*(eta+beta*nu)'*hat3(w)*(hat3(w)-alpha*rot_log(eye(3),RRef))*(zeta+alpha*nu);
% Compute <m2*JYh+m3*JYv,D_JYh_JXv>
g_D_JYh_Jxv = m3/2*(eta+beta*nu)'*(kv*hat3(w)+beta*rot_log(eye(3),RRef)+kd*RRef'*rot_log(RRef,R))*(zeta+alpha*nu)...
    -( m2*(zeta+alpha*nu)+m3*(eta+beta*nu) )'*kd*rot3_logDiffMat(RRef,R)*(zeta+alpha*nu)...
    +( m2*(zeta+alpha*nu)+m3*(eta+beta*nu) )'*(kd*rot3_logDiffMat(R,RRef)-beta*rot3_logDiffMat(eye(3),RRef))*nu...
    -( m2*(zeta+alpha*nu)+m3*(eta+beta*nu) )'*kv*(eta+beta*nu);

% Sum all parts to compute the metric
g_DJYJX_JY = g_D_Y_X_SO3 + g_DJYh_JXh + g_R_U_JYh_JXh + g_R_JYh_JXh_U...
    + g_D_JYh_Jxv;
end

