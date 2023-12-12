function [D, output] = rotBundle_metric_nonNatural_kozul(U, X, Y, Z, m, varargin)
% Implements the non-natural metric resulting from the Kozul Formula for
% TSO(3). IE: D = 2*g_bar(Del_{X}Y, Z), where g_bar is the metric defined
%   g_bar(X^h, Y^h) = m1*g(X,Y)
%   g_bar(X^h, Y^v) = m2*g(X,Y)
%   g_bar(X^v, Y^v) = m3*g(X,Y)
% where g is the metric on the SO(3)
% INPUTS:
%   U(R) := the tangent vector at R which creates the curve on TSO(n) in
%       the form of [R; U]
%   X(R,U), Y(R,U), Z(R,U) := Vector fields on TSO(3) given as matrices 
%       (IE [6 x 3])
%   m := an array of 3 values for the non natural metric on TSO(3)
% OUTPUTS:
%   D(R) := Result of g_bar(Del_{X}Y, Z) on TSO(3)

flagM3Terms=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'extras'
            flagM3Terms=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end


%% Standard Variables
if (length(m) ~= 3)
    error('m must be an array of length 3')
end
m1 = m(1); m2 = m(2); m3 = m(3); % These should produce a pos. def. symm. matrix M

%% Split X, Y, Z into vertical and horizontal components
PU = @(R) [R; U(R)];
Xh = @(R) rotBundle_extractHoriz(PU(R), X(R, U(R)));
Xv = @(R) rotBundle_extractVert(PU(R), X(R, U(R)));
Yh = @(R) rotBundle_extractHoriz(PU(R), Y(R, U(R)));
Yv = @(R) rotBundle_extractVert(PU(R), Y(R, U(R)));
Zh = @(R) rotBundle_extractHoriz(PU(R), Z(R, U(R)));
Zv = @(R) rotBundle_extractVert(PU(R), Z(R, U(R)));

%% Compute all combinations of the vertical and horizontal components
% Notation of variable names. gbar_Del_X_Y_Z means take covariant of Y wrt
% X then take inner product of that with Z.
gbar_Del_Xh_Yh_Zh = @(R) m1*rot_metric(R, rot_covar(Xh, Yh), Zh) + ...
    m2*rot_metric(R, rot_curvature(R, U, Xh, Yh), Zh);
gbar_Del_Xh_Yh_Zv = @(R) m2*rot_metric(R, rot_covar(Xh, Yh), Zv)- ...
    m3/2*rot_metric(R, Zv, rot_curvature(R, Xh, Yh, U));
gbar_Del_Xh_Yv_Zh = @(R) m2*rot_metric(R, Zh, rot_covar(Xh, Yv)) + ...
    m3/2*rot_metric(R, rot_curvature(R, U, Yv, Xh), Zh);
gbar_Del_Xh_Yv_Zv = @(R) m3*rot_metric(R, Zv, rot_covar(Xh, Yv));
gbar_Del_Xv_Yh_Zh = @(R) m3/2*rot_metric(R, rot_curvature(R, U, Xv, Yh), Zh);
gbar_Del_Xv_Yh_Zv = @(R) 0;
gbar_Del_Xv_Yv_Zh = @(R) 0;
gbar_Del_Xv_Yv_Zv = @(R) 0;

%% Resulting metric is the sum of all the components above
D = @(R) gbar_Del_Xh_Yh_Zh(R) + gbar_Del_Xh_Yh_Zv(R) + ...
    gbar_Del_Xh_Yv_Zh(R) + gbar_Del_Xh_Yv_Zv(R) + ...
    gbar_Del_Xv_Yh_Zh(R) + gbar_Del_Xv_Yh_Zv(R) + ...
    gbar_Del_Xv_Yv_Zh(R) + gbar_Del_Xv_Yv_Zv(R);

if flagM3Terms
    output.gbar_Del_Xh_Yh_Zh= gbar_Del_Xh_Yh_Zh;
    output.gbar_Del_Xh_Yh_Zv= gbar_Del_Xh_Yh_Zv;
    output.gbar_Del_Xh_Yv_Zh = gbar_Del_Xh_Yv_Zh;
    output.gbar_Del_Xh_Yv_Zv = gbar_Del_Xh_Yv_Zv;
    output.gbar_Del_Xv_Yh_Zh = gbar_Del_Xv_Yh_Zh;    
end

