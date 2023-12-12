function [M_ALL] = rotRef_contractionMat2(beta, kd, kv, kref, M_nonnatural, varargin)
% LAST EDIT: OCT 4, 2020 (change kp to kref to match CDC 2020 paper)
% Computes the matrix form of the M matrix (on TSO(3)xSO(3)) such that 
% [zeta;eta;nu]'*M*[zeta;eta;nu] == g(\nabla_{Vx}X + lambda*Vx,Vx)_m, where
% zeta and eta are the R^3 representation of the horizontal and vertical
% tangent vectors of Vx on TSO(3) and nu is the R^3 vector on SO(3)
% ASSUMES: Vx = [R*hat3(zeta);R*hat3(eta);RRef*hat3(nu)];
% ASSUMES: X = [U;kd*rot_log(R,RRef)-kv*U;kref*log(RRef,eye(3))];
% NOTE: When calling function handle, w should be a [3x1] vector
% INPUTS:
%   beta := scalar for the exponential convergence rate
%   kd := scalar positive position gain
%   kv := scalar positive velocity gain
%   kref := scalar positive reference traj. gain
%   M_nonnatural := nonnatural metric gain matrix on TSO(3)xSO(3) in the
%       form of [m1 m2 m6;m2 m3 m5;m6 m5 m4]
% OUTPUTS:
%   M_ALL(R,w,RRef) := M matrix such that 
%       [zeta;eta;nu]'*M*[zeta;eta;nu] == g(\nabla_{Vx}X + beta*Vx,Vx)_m

flagOutputSymmetric = false;

%optional parameters
ivarargin=1;
len_varargin = length(varargin);
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'sym'
            flagOutputSymmetric = true;
        otherwise
            % ignore the entry and increase index until next entry is
            % a string
            while ivarargin+1 <= len_varargin && ~ischar(varargin{ivarargin+1})
                ivarargin=ivarargin+1;
            end 
    end
    ivarargin=ivarargin+1;
end

%% Standard Variables
% Define gains for metric on TSO(3)xSO3
% Transform the metric from the "nonnatural" coordinates into the "natural"
% coordinates so we can compute the contraction matrix using the natural
% covar on TSO3xSO3 using Schur's complement
% [J,M_natural] = rotRef_SchurComplement2(M_nonnatural);
m1 = M_nonnatural(1,1); m2 = M_nonnatural(1,2); m3 = M_nonnatural(2,2);
m4 = M_nonnatural(3,3); m5 = M_nonnatural(2,3); m6 = M_nonnatural(1,3);
m1_prime = m1-m6^2/m4; m2_prime = m2-m5*m6/m4; m3_prime = m3-m5^2/m4;

%% Compute the components of the M matrix on TSO(3) (THIS SHOULD BE DIFFERENT DUE TO log(R,RREF) term!!!
M11 = @(R,w,RRef) m2_prime/4*hat3(w)^2 ...
    -m2*kd*rot3_logDiffMat(RRef,R);

M21 = @(R,w,RRef) -m2_prime/4*hat3(w)...
    + m3_prime/8*hat3(w)^2 ...
    + m3_prime*kv/4*hat3(w)...
    + m3_prime*kd/4*rot_log(eye(3),RRef'*R)...
    -m3*kd/2*rot3_logDiffMat(RRef,R)...
    +(m1-m2*kv)/2*eye(3);
M22 = @(R,w,RRef) (m2-m3*kv)*eye(3);

%% Compute the cross and SO3 components
M31 = @(R,w,RRef) m2*kd/2*rot3_logDiffMat(RRef,R)...
    -m5*kd/2*rot3_logDiffMat(RRef,R)...
    +m6/4*(-kref*rot_log(eye(3),RRef)+m6/m4*hat3(w)+m5/m4*kd*rot_log(eye(3),R'*RRef)-m5/m4*kv*hat3(w))...
    -m6/2*kref*rot3_logDiffMat(RRef,eye(3));
M32 = @(R,w,RRef) m3*kd/2*rot3_logDiffMat(RRef,R)...
    +(m6-m5*kv)/2*eye(3)...
    +m5/4*(-kref*rot_log(eye(3),RRef)+m6/m4*hat3(w)+m5/m4*kd*rot_log(eye(3),R'*RRef)-m5/m4*kv*hat3(w))...
    -m5/2*kref*rot3_logDiffMat(RRef,eye(3));
M33 = @(R,w,RRef) -m4*kref*rot3_logDiffMat(eye(3),RRef)...
    +m5*kd*rot3_logDiffMat(R,RRef);

%% Construct the complete matrix by adding new terms from the reference manifold
M_ALL = @(R,w,RRef) [M11(R,w,RRef) M21(R,w,RRef)' M31(R,w,RRef)';...
    M21(R,w,RRef) M22(R,w,RRef) M32(R,w,RRef)';...
    M31(R,w,RRef) M32(R,w,RRef) M33(R,w,RRef)] + kron(beta*M_nonnatural,eye(3));
    
if flagOutputSymmetric
    M_ALL = @(R,w,RRef) (M_ALL(R,w,RRef) + M_ALL(R,w,RRef)')/2;
end
end

