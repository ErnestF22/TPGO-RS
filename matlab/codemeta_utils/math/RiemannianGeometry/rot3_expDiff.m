%function BexpA=rot3_expDiff(R,A,B)
%Computes the differential (dExp_R|A) applied to B. In other words, given a line l(t)=A+t*B in
%T_R SO(3) gives the tangent vector in T_{Exp_R(A)} SO(3) for the curve
%Exp_R(l(t))
function DA=rot3_expDiff(R,A,B,varargin)
%compute result by going to local coordinates, applying the matrix version,
%and then go back to tangent vectors
BVec=rot_vee(R,B);
[D,R2]=rot3_expDiffMat(R,A,varargin{:});
DAVec=D*BVec;
DA=rot_hat(R2,DAVec);
