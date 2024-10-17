%function vA=rot3_expDiffInv(R,A,BexpA)
%Computes the inverse of the differential (dExp_R|A) applied to BexpA
%(node: BexpA is a tangent vector at Exp_R(A)). Given the tangent to the
%curve passing through Exp_R(A) with direction BexpA, this function gives
%the vector B in T_SO(3)R such that the curve Exp_R(A+t*B) has the same
%tangent as the original curve
function B=rot3_expDiffInv(R,A,BexpA,varargin)
flagAisRot=false;   %option to say that A, in fact, is Exp_R(A)

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'rot'
            flagAisRot=true;
    end
    ivarargin=ivarargin+1;
end

if flagAisRot
    R2=A;
    A=rot_log(R,R2);
else
    R2=rot_exp(R,A);
end
R2=R'*R2;

%pull back to the identity
A=R'*A;
BexpA=R'*BexpA;

%compute and then push forward
B=R*rot_hat(eye(3),rot3_expDiffSeries(A)\rot_vee(eye(3),R2'*BexpA));
