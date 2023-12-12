function [M_ALL] = rotRefOpt_contractionMat_kdot(varargin)
% LAST EDIT: OCT 28, 2020 by Bee Vang
% THIS IS A COPY OF rotRef_contractionMat2.m on Oct. 4, 2020.
% Computes the derivative wrt [kd,kv,kref] of the matrix form of the M matrix (on TSO(3)xSO(3)) such that 
% [zeta;eta;nu]'*M*[zeta;eta;nu] == g(\nabla_{Vx}X + lambda*Vx,Vx)_m, where
% zeta and eta are the R^3 representation of the horizontal and vertical
% tangent vectors of Vx on TSO(3) and nu is the R^3 vector on SO(3)
% ASSUMES: Vx = [R*hat3(zeta);R*hat3(eta);RRef*hat3(nu)];
% ASSUMES: X = [U;kd*rot_log(R,RRef)-kv*U;kref*log(RRef,eye(3))];
% NOTE: When calling function handle, w should be a [3x1] vector
% NOTE: In our optimal control approach, we use this function to compute
%   the quadratic part of the contraction metric (all known) which doesn't
%   involve the rate of change of the gains. This function handle should
%   only be created once and used multiple times to recompute the quadratic
%   part at every (R,w,Rref,kd,kv,kref)
% INPUTS:
%   M_nonnatural := nonnatural metric gain matrix on TSO(3)xSO(3) in the
%       form of [m1 m2 m6;m2 m3 m5;m6 m5 m4]
% OUTPUTS:
%   M_ALL(R,w,RRef,gains_dot,M_nn) := Contraction matrix with arguements
%       R := current rotation [3x3]
%       w := current angular velocity [3x1]
%       RRef := reference rotation [3x3]
%       gains_dot := der vector of positive system gains [kd;kv;kref] [3x1]
%       M_nn := non-natural metric on TSO(3)xSO(3) [m1,m2,m6;m2,m3,m5;m6,m5,m4] [3x3]

%% Load Optional Parameters
flagOutputSymmetric = false;
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

%% Compute the components of the M matrix on TSO(3) (THIS SHOULD BE DIFFERENT DUE TO log(R,RREF) term!!!
M11 = @(R,w,RRef,gains_dot,M_nn) -M_nn(1,2)*gains_dot(1)*rot3_logDiffMat(RRef,R);
M21 = @(R,w,RRef,gains_dot,M_nn)(M_nn(2,2)-M_nn(2,3)^2/M_nn(3,3))*gains_dot(2)/4*hat3(w)...
    + (M_nn(2,2)-M_nn(2,3)^2/M_nn(3,3))*gains_dot(1)/4*rot_log(eye(3),RRef'*R)...
    -M_nn(2,2)*gains_dot(1)/2*rot3_logDiffMat(RRef,R)...
    -M_nn(1,2)*gains_dot(2)/2*eye(3);
M22 = @(R,w,RRef,gains_dot,M_nn) -M_nn(2,2)*gains_dot(2)*eye(3);
%% Compute the cross and SO3 components
M31 = @(R,w,RRef,gains_dot,M_nn) M_nn(1,2)*gains_dot(1)/2*rot3_logDiffMat(RRef,R)...
    -M_nn(2,3)*gains_dot(1)/2*rot3_logDiffMat(RRef,R)...
    +M_nn(1,3)/4*(-gains_dot(3)*rot_log(eye(3),RRef)+M_nn(2,3)/M_nn(3,3)*gains_dot(1)*rot_log(eye(3),R'*RRef)-M_nn(2,3)/M_nn(3,3)*gains_dot(2)*hat3(w))...
    -M_nn(1,3)/2*gains_dot(3)*rot3_logDiffMat(RRef,eye(3));
M32 = @(R,w,RRef,gains_dot,M_nn) M_nn(2,2)*gains_dot(1)/2*rot3_logDiffMat(RRef,R)...
    -M_nn(2,3)*gains_dot(2)/2*eye(3)...
    +M_nn(2,3)/4*(-gains_dot(3)*rot_log(eye(3),RRef)+M_nn(2,3)/M_nn(3,3)*gains_dot(1)*rot_log(eye(3),R'*RRef)-M_nn(2,3)/M_nn(3,3)*gains_dot(2)*hat3(w))...
    -M_nn(2,3)/2*gains_dot(3)*rot3_logDiffMat(RRef,eye(3));
M33 = @(R,w,RRef,gains_dot,M_nn) -M_nn(3,3)*gains_dot(3)*rot3_logDiffMat(eye(3),RRef)...
    +M_nn(2,3)*gains_dot(1)*rot3_logDiffMat(R,RRef);

%% Construct the complete matrix by adding new terms from the reference manifold
M_ALL = @(R,w,RRef,gains_dot,M_nn) [M11(R,w,RRef,gains_dot,M_nn) M21(R,w,RRef,gains_dot,M_nn)' M31(R,w,RRef,gains_dot,M_nn)';...
    M21(R,w,RRef,gains_dot,M_nn) M22(R,w,RRef,gains_dot,M_nn) M32(R,w,RRef,gains_dot,M_nn)';...
    M31(R,w,RRef,gains_dot,M_nn) M32(R,w,RRef,gains_dot,M_nn) M33(R,w,RRef,gains_dot,M_nn)];

if flagOutputSymmetric
    M_ALL = @(R,w,RRef,gains_dot,M_nn) (M_ALL(R,w,RRef,gains_dot,M_nn) + M_ALL(R,w,RRef,gains_dot,M_nn)')/2;
end
end

