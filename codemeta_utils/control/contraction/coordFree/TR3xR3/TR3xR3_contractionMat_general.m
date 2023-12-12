function [M_contract] = TR3xR3_contractionMat_general(hessf_R_Rd, hessf_Rd_0, kd,kv,kp,beta,M_nn)
% Last Edited: Oct. 27th 2020 by Bee Vang
% Compute the matrix form of the M matrix (on TR^3xR^3) such that
% [zeta;eta;nu]'*M*[zeta;eta;nu] = <D_Y_X+lambda*Y,Y> where zeta and eta
% are the R^3 representation of the horizontal and vertical tangent vectors
% of Y on TR^3 and nu is the R^3 vector on R^3. Requires the hessian of the
% chosen cost functions \rho_x and \rho_{x_d}
% INPUTS:
%   hessf_R_Rd(A,CC_R_Rd) := function handle for hessian of cost from current to
%      reference position
%   hessf_Rd_0(C,~) := function handle for hessian of cost from reference
%      position to the identity
%   kd,kv,kp := scalar [position;velocity;reference position] gains
%   beta := positive scalar min exponential rate
%   M_nn := A [3x3] pos. def. matrix representing the metric gains.
% OUTPUTS:
%   M_contract := M matrix such that 
%       [zeta;eta;nu]'*M*[zeta;eta;nu] = <D_Y_X+lambda*Y,Y>
%       NOTE: This should be a function of (A,B,C) but it doesnt matter in
%       R^6
% Extract the gains
m1 = M_nn(1,1); m2 = M_nn(1,2); m3 = M_nn(2,2);
m4 = M_nn(3,3); m5 = M_nn(2,3); m6 = M_nn(1,3);
% Define the hessian
% hessf = @(x,xd) (sqrt(1+norm(x-xd)^2/2)*eye(3)-(x-xd)*(x-xd)'/(2*sqrt(1+norm(x-xd)^2/2)))/(sqrt(1+norm(x-xd)^2/2))^2;
% Define the components
M11 = @(A,B,C) -2*m2*kd*hessf_R_Rd(A,C);
M21 = @(A,B,C) -m3*kd*hessf_R_Rd(A,C)+(m1-m2*kv)*eye(3);
M22 = @(A,B,C) (2*m2-2*m3*kv)*eye(3);
M31 = @(A,B,C) -m5*kd*hessf_R_Rd(A,C)...
    +m2*kd*hessf_R_Rd(A,C)-m6*kp*hessf_Rd_0(C,zeros(3,1));
M32 = @(A,B,C) m3*kd*hessf_R_Rd(A,C)...
    -m5*kp*hessf_Rd_0(C,zeros(3,1))+(m6-m5*kv)*eye(3);
M33 = @(A,B,C) 2*m5*kd*hessf_R_Rd(A,C)-2*m4*kp*hessf_Rd_0(C,zeros(3,1));

M_contract = @(A,B,C) [M11(A,B,C) M21(A,B,C)' M31(A,B,C)'; ...
    M21(A,B,C) M22(A,B,C) M32(A,B,C); ...
    M31(A,B,C) M32(A,B,C) M33(A,B,C)] + kron(beta*M_nn,eye(3));
end

