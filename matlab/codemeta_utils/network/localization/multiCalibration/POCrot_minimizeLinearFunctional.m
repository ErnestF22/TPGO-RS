%Minimize linear functional on O(n)
%function rot_minimizeLinearFunctional(A,b,R0,varargin)
%Minimizes a functional of the form norm(A*R(:)-b,2)
function [R,output]=rot_minimizeLinearFunctional(A,b,R0,varargin)
tolNorm=1e-14;
maxIt=100;
flagCollectOutput=false;
stepSize=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'tolnorm'
            ivarargin=ivarargin+1;
            tolNorm=varargin{ivarargin};
        case 'stepsize'
            ivarargin=ivarargin+1;
            stepSize=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if nargout>1
    flagCollectOutput=true;
end

d=size(R0,1);
EHat=reshape(rot_tangentBasis(eye(d)),d^2,[]);
R=R0;
cprev=Inf;
for it=1:maxIt
    f=A*R(:)-b;
    c=f'*f;
    if c>cprev+1e-9
        warning('Cost went up at iteration %d. Stopping here.',it)
        break;
    end
    v=(A*kron(eye(d),R)*EHat)\f;
    R=rot_exp(R,rot_hat(R,stepSize*v));
    if norm(v)<tolNorm
        break;
    end
    cprev=c;
end

if flagCollectOutput
    output.c=c;
    output.it=it;
    output.v=v;
    output.normv=norm(v);
end
