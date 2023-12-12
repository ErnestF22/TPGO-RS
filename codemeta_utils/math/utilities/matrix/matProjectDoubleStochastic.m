%Project a matrix on the set of double stochastic matrices
%function [X,output]=projectDoubleStochastic(Z)
%Projects Z to the space of double stochastic matrices under various norms
%Inputs
%   Z       matrix to project
%Optional inputs
%   {'fro','l1','entropy',normalized'}  string indicating what distance to minimize
%       'fro'           (default) Frobenious norm of difference
%       'l1'            L-1 norm (ratio cut)
%       'normalized'    normalize with degree (normalized cut)
%       'entropy'       relative entropy (iterated 'normalized')
%       
%   'tol',tol   Tolerance for entropy method
%
% For proofs of the formulas, see 
%   R. Zass, A. Shashua. 
%   "Doubly stochastic normalization for spectral clustering."
%   Advances in Neural Information Processing Systems (NIPS), 2006.
function [X,output]=matProjectDoubleStochastic(Z,varargin)
method='fro';
tol=1e-12;
maxIt=500;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case {'fro','l1','normalized','entropy'}
            method=lower(varargin{ivarargin});
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch lower(method)
    case 'fro'
        n=size(Z,1);
        U=1/n*ones(n);
        P=eye(n)-U;
        X=P*Z*P+U;
        output.P=P;
        output.U=U;
    case 'l1'
        D=eye(size(Z))-diag(diag(Z));
        X=Z+D;
        output.D=D;
    case 'normalized'
        [X,DSqrtInv]=degreeNormalize(Z);
        output.DSqrtInv=DSqrtInv;
    case 'entropy'
        X=Z;
        DSqrtInv=eye(size(X));
        for it=1:maxIt
            [X,DSqrtInvCur,e]=degreeNormalize(X);
            DSqrtInv=DSqrtInv*DSqrtInvCur;
            if abs(e)<tol
                break
            end
        end
        output.DSqrtInv=DSqrtInv;
    otherwise
        error('Not implemented')
end

function [X,DSqrtInv,e]=degreeNormalize(X)
d=sum(X,2);
e=norm(d-1,'inf');

dSqrt=sqrt(d);
DSqrtInv=diag(1./dSqrt);
X=DSqrtInv*X*DSqrtInv;
