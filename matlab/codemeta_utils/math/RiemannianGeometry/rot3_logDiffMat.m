%Compute matrix representation of the differential of the logarithm in SO(3)
%NOTE: This function assumes that only A is changing (see rot3_logDiff.m)
%function D=rot3_logDiff(R,A,varargin)
%Inputs
%   R   indicates in which tangent space the differential should be
%       computed
%   A   indicates to which point the space the differential should be
%       computed (by default, this needs to be rotation, but see also the
%       optional arguments)
%Optional arguments
%   'rot'           A is a rotation
%   'tangent'       A is a [3x3] matrix containing a tangent vector at R
%   'tangentVec'    A is a [3x1] vector containing a tangent vector in
%                   coordinates
function D=rot3_logDiffMat(R,A,varargin)
tolTheta=1e-12;
flagAisRot=true;            %option to say that A, in fact, is Exp_R(A)
flagAisTangentVec=false;    %option to say that A is in vector coordinate representation

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'rot'
            flagAisRot=true;
        case 'tangent'
            flagAisRot=false;
            flagAisTangentVec=false;
        case 'tangentVec'
            flagAisRot=false;
            flagAisTangentVec=true;
    end
    ivarargin=ivarargin+1;
end

if flagAisRot
    R2=A;
    A=rot_log(R,R2);
end

if ~flagAisTangentVec
    %pull back to the identity
    A=R'*A;
    vecA=vee3(A);
%     vecA = rot_vee(R, A);
else
    vecA=A;
end
[u,theta]=cnormalize(vecA);

if theta<tolTheta
    D=eye(3);
else
    uut=u*u';

    D=uut+theta*(hat3(u)-cot(theta/2)*hat3(u)^2)/2;
end