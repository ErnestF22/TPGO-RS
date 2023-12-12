function [Ft, dFt, appder] = TR3xR3_metric_compatibility_nonnatural(A,B,C,t,X,Y,Z,M_nn,varargin)
% Check if d/dt<X,Y> = <D_{Z}X,Y> + <X,D_{Z}Y> evaluated at (A,B,C,t).
% Ideally Z=Y.
% INPUTS:
%   A(t) := The evaluation position on TR^3
%   B(A,t) := The evaluation velocity on TR^3
%   C(t) := The evaluation position on R^3
%   t := Time of evaluation
%   X(A,B,C), Y(A,B,C), Z(A,B,C) := Tangent Vector as [9x1] vector.
%       NOTE: Z=d/dt Z_curve
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.
% OUTPUTS:
%   Ft := numerical evaluation of <X,Y>
%   dFt := analytical evaluation of <D_{Z}X, Y> + <X, D_{Z}Y>
%   appder := numerical approximation of d/dt <X,Y>

% Define useful parameters
Z_curve = @(t) [A(t);B(A(t),t);C(t)];
B_A = @(A) B(A,t);

% Optional Parameters
flagContractionMatrix_2norm = false;
flagContractionMatrix_general = false;
kd = 0;
kv = 0;
kp = 0;
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'contractionmatrix_2norm'
            % flag to test <D_Z_X,Y> using the contraction matrix, requires
            % kd and kv gains
            flagContractionMatrix_2norm = true;
            ivarargin=ivarargin+1;
            kd = varargin{ivarargin};
            ivarargin=ivarargin+1;
            kv = varargin{ivarargin};
        case 'contractionmatrix_general'
            flagContractionMatrix_general = true;
            ivarargin=ivarargin+1;
            gradf_A_C = varargin{ivarargin};
            ivarargin=ivarargin+1;
            gradf_C_0 = varargin{ivarargin};
            ivarargin=ivarargin+1;
            kd = varargin{ivarargin};
            ivarargin=ivarargin+1;
            kv = varargin{ivarargin};
            ivarargin=ivarargin+1;
            kp = varargin{ivarargin};
        otherwise    
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


% Define the metric as a function of t
Xt = @(t) X(A(t),B(A(t),t),C(t));
Yt = @(t) Y(A(t),B(A(t),t),C(t));
Left = @(t) TR3xR3_metric_nonnatural(Xt(t),Yt(t),M_nn);

% Compute the right side <D_Z_X,Y>+<X,D_Z_Y>
% <D_Z_X,Y>
if flagContractionMatrix_2norm
    contractionMat = TR3xR3_contractionMat(kd,kv,M_nn);
    Y_numeric = Y(A(t),B(A(t),t),C(t));
    zeta = Y_numeric(1:3); eta = Y_numeric(4:6); nu = Y_numeric(7:9);
    Yvec = [zeta;eta;nu];
    g_DZX_Y = @(t) Yvec'*contractionMat*Yvec;
elseif flagContractionMatrix_general
    contractionMat = TR3xR3_contractionMat_general(gradf_A_C,...
        gradf_C_0, kd, kv, kp, 0, M_nn);
    Y_numeric = Y(A(t),B(A(t),t),C(t));
    zeta = Y_numeric(1:3); eta = Y_numeric(4:6); nu = Y_numeric(7:9);
    Yvec = [zeta;eta;nu];
    g_DZX_Y = @(t) Yvec'*contractionMat(A(t),B(A(t),t),C(t))*Yvec/2;
    %NOTE: we divide by 2 above because the original contraction paper 
    %"On Contraction Analysis for Nonlinear Systems" (Lohmiller & Slotine)
    %multiply both sides by 2
else
    D_Z_X = TR3xR3_covar_directMethod(A,B,C,t,Z,X,M_nn);
    g_DZX_Y = @(t) TR3xR3_metric_nonnatural(D_Z_X,Yt(t),M_nn);
    % check if direct result matches the coord change method
    [D_Z_X_nonnatural,D_Z_X_natural] = TR3xR3_covar_coordChange(A,B,C,t,Z,X,M_nn);
    if any(abs(D_Z_X_nonnatural-D_Z_X)>1e-6)
        error('error computing covar using coordinate change method')
    end
end

% <D_Z_Y,X>
D_Z_Y = TR3xR3_covar_directMethod(A,B,C,t,Z,Y,M_nn);
g_DZY_Y = @(t) TR3xR3_metric_nonnatural(D_Z_Y,Xt(t),M_nn);
% Combine terms to compute full derivative
g_nonnatural_total_der_fun = @(t) g_DZX_Y(t) + g_DZY_Y(t);
% Perform checek
[Ft, dFt, appder] = funCheckDer(Left, g_nonnatural_total_der_fun, t, 'nodisplay');
end