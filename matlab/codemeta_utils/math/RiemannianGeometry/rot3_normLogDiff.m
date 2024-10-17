function D=rot3_normLogDiff(R,A,varargin)
flagAisRot=true;   %option to say that A, in fact, is Exp_R(A)

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'rot'
            flagAisRot=true;
        case 'tangent'
            flagAisRot=false;
    end
    ivarargin=ivarargin+1;
end

if flagAisRot
    R2=A;
    A=rot_log(R,R2);
end

%pull back to the identity
A=R'*A;

[u,theta]=cnormalize(vee(A));
uut=u*u';

D=(hat(u)-cot(theta/2)*hat(u)^2)/2;
