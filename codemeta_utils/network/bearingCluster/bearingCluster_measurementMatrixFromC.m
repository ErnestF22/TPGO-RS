function M=bearingCluster_measurementMatrixFromC(C,u,varargin)
methodNormalization='orthogonal';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'methodnormalization'
            ivarargin=ivarargin+1;
            methodNormalization=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

U=bearingCluster_augmentedBearingMatrix(u);

d=size(U,1)/size(U,2);
Cd=kron(C,eye(d));
M=Cd*U;
switch lower(methodNormalization)
    case 'orthogonal'
        [U,S,V]=svd(Cd*Cd');
        Dinv=U*sqrt(S)*U';
    case 'cholesky'
        Dinv=chol(Cd*Cd');
    case 'none'
        Dinv=eye(size(Cd,1));
    otherwise
        error('Normalization method not recognized')
end

M=Dinv\M;
