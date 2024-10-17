%function D=rot3_expDiffInvMat(R,A,varargin)
%Gives the 3x3 matrix representation of rot3_expDiffInv(R,A), where the
%basis for the tangent spaces are given by rot_tangentBasis
function D=rot3_expDiffInvMat(R,A,varargin)
flagAisRot=true;   %option to say that A, in fact, is Exp_R(A)
flagMethod='closedForm';    %method to use: can be 'closedForm' or 'series'

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'rot'
            flagAisRot=true;
        case 'tangent'
            flagAisRot=false;
        case 'method'
            ivarargin=ivarargin+1;
            flagMethod=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid'])
    end
    ivarargin=ivarargin+1;
end

if flagAisRot
    R2=A;
    A=rot_log(R,R2);
else
    R2=rot_exp(R,A);
end

%pull back to the identity
A=R'*A;

switch lower(flagMethod)
    case 'inverseseries'
        D=inv(rot3_expDiffMat(eye(3),A,'series'));
    case 'inverseclosedform'
        D=inv(rot3_expDiffMat(eye(3),A,'closedform'));
    case 'closedform'
        [u,theta]=cnormalize(vee3(A));
        uut=u*u';

        D=(uut)'+theta/2*cot(theta/2)*(eye(3)-uut)+theta/2*hat3(u);
    otherwise
        error('Method not valid')
end
