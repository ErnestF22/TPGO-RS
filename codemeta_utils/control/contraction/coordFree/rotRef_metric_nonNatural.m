function [D] = rotRef_metric_nonNatural(Z,X,Y,M)
% Compute the non-natural metric on TSO(3)xSO(3) with scaling factors m_{i}
% INPUTS:
%   Z := Point of evaluation (p,u,p') given as a matrix (IE [9x3]) on
%       TSO(3)xSO(3)
%   X,Y := Vector fields on TSO(3)xSO(3)
%   M := A matrix of scaling factors (should be symmetric), [3x3] matrix
% OUTPUTS:
%   D := result of <X,Y> on TSO(3)xSO(3)

% Define Parameters
% Scaling factors for metric on TSO(3)
m1 = M(2,2); m2 = M(2,3); m3 = M(3,3);
% Scaling factors for metric on SO(3) (the reference trajectory)
m4 = M(1,1); % g(X_SO3,Y_SO3)
% Scaling factors for cross metric on TSO(3)xSO(3)
m5 = M(1,2); % g(X^h,Y_SO3)
m6 = M(1,3); % g(X^v,Y_SO3)

% Compute metric on TSO(3)
Z_TSO3 = Z(4:9,:);X_TSO3 = X(4:9,:);Y_TSO3 = Y(4:9,:);
M_TSO3 = [m1;m2;m3];
D_TSO3 = rotBundle_metric_nonNatural(Z_TSO3, X_TSO3, Y_TSO3, M_TSO3);

% Compute metric on SO(3) (the reference trajectory)
R_SO3 = Z(1:3,:);X_SO3 = X(1:3,:);Y_SO3 = Y(1:3,:);
D_SO3 = m4*rot_metric(R_SO3,X_SO3,Y_SO3);
% Left translate from RRef to R
R_TSO3 = Z(4:6,:);
X_SO3_translate = R_TSO3*R_SO3'*X_SO3;
Y_SO3_translate = R_TSO3*R_SO3'*Y_SO3;
% Compute cross (manifold) terms (computed at R_TSO3)
X_h = X(4:6,:); X_v = X(7:9,:);
Y_h = Y(4:6,:); Y_v = Y(7:9,:);
D_Xh_YSO3 = m5*rot_metric(R_TSO3,X_h,Y_SO3_translate);
D_Xv_YSO3 = m6*rot_metric(R_TSO3,X_v,Y_SO3_translate);
D_XSO3_Yh = m5*rot_metric(R_TSO3,X_SO3_translate,Y_h);
D_XSO3_Yv = m6*rot_metric(R_TSO3,X_SO3_translate,Y_v);

% The resulting metric is a combination of the metric on TSO3 and SO3
D = D_TSO3+D_SO3 + D_Xh_YSO3 + D_Xv_YSO3 + D_XSO3_Yh + D_XSO3_Yv;
end

