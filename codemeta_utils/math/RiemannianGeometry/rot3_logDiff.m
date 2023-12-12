function [D] = rot3_logDiff(R,A,dA)
% Computes the matrix representation of the differential of the logarthm in SO(3)
% NOTE: this takes into account that only A is changing
% NOTE: Useful relationship log_{R}(A) = -RA'*log_{A}(R)
% INPUTS:
%   R  := indicates in which tangent space the differential should be
%       computed
%   dR := the tangent vector associated with the change along R
%   A  := indicates to which point the differential should be
%       computed (by default, this needs to be rotation, 
%       TODO: optional arguments)
%   dA := the tangent vector associated with the change along A

%%% If we want to consider that R is also changing, we can add these two
%%% terms. For now, we're leaving them out since they can be easily
%%% implented outside of this function and this implementation allows for
%%% greater flexibility.
% % Consider how the log changes as R changes, resulting in two terms
% % following the relationship above.
% % First term is how R is changing outside of the log
% d_log_A_R_Rdot = -dR*A'*rot_log(A,R);
% % Second term is how R is changing inside the log
% d_log_A_R_DlogR = -rot_hat(R*A', rot3_logDiffMat(A,R)*rot_vee(R,dR));
%%%
d_log_A_R_Rdot = 0;
d_log_A_R_DlogR = 0;
d_log_R_A_DlogA = rot_hat(R, rot3_logDiffMat(R,A)*rot_vee(A,dA));


% Result is the combination of the three terms
D = d_log_A_R_Rdot + d_log_A_R_DlogR + d_log_R_A_DlogA;


end

