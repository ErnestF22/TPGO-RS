function [X,output] = sfm_rawRotationsLowRankALM_SolverL1(M,flag,varargin)

% this program solves
% min ||M-X||_1 + lam ||X||_*
% or
% min ||M-X||_1 s.t. X is PSD
% st: X is symmetric, each sub-block of X is in SO(3), diagonal are identities
n = size(M,1)/3;
lambda = sqrt(3*n); % the regularization weight
lambdaFactor=0.5;
innerProj=true;
groupSparse=false;
M0=M;
tol=1e-4;
maxIter=4000;
flagShowStats=false;
flagCollectStats=nargout>1;
methodLowRankPrior='nuclear';

output=[];

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        case 'lambdafactor'
            ivarargin=ivarargin+1;
            lambdaFactor=varargin{ivarargin};
        case 'innerproj'
            ivarargin=ivarargin+1;
            innerProj=varargin{ivarargin};
        case 'lowrankprior'
            ivarargin=ivarargin+1;
            methodLowRankPrior=lower(varargin{ivarargin});
        case 'groupsparse'
            ivarargin=ivarargin+1;
            groupSparse=varargin{ivarargin};
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        case 'maxiter'
            ivarargin=ivarargin+1;
            maxIter=varargin{ivarargin};
        case 'showstats'
            flagShowStats=true;
        case 'collectstats'
            flagCollectStats=true;
        case 'init'
            ivarargin=ivarargin+1;
            M0=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

lambda=lambdaFactor*lambda;

%% initialization
S = false(size(M));
for i = 1:n-1
    for j = i+1:n
        if flag(i,j)
            S(3*i-2:3*i,3*j-2:3*j) = true(3,3);
        end
    end
end

X = M0;
% Z = M;
% E = zeros(size(M));

Y1 = zeros(size(M)); % the dual variable
Y2 = Y1;
mu1 = .25/mean(abs(M(:)));
mu2 = mu1;

if flagCollectStats
    output.maxIter=maxIter;
    output.primRes=NaN(maxIter,2);
    output.dualRes=NaN(maxIter,2);
    output.mu=NaN(maxIter,2);
    output.lambda=lambda;
end

for iter = 1:maxIter
    
    % update E
    E = (M-X+Y1/mu1).*S;
    if groupSparse
        for i = 1:n-1
            for j = i+1:n
                if flag(i,j)
                    E(3*i-2:3*i,3*j-2:3*j) = groupthresholdL2(E(3*i-2:3*i,3*j-2:3*j),1/mu1);
                end
            end
        end
    else
        E = sign(E).*max(abs(E)-1/mu1,0);
    end
    
    % update Z
    switch methodLowRankPrior
        case 'nuclear'
            [U,W,V] = svd(X + Y2/mu2);
            w = diag(W) - lambda/mu2;
            ind = find(w>0);
            Z = U(:,ind)*diag(w(ind))*V(:,ind)';
        case 'sdp'
            try
                [U,W]=eig(X+Y2/mu2);
            catch
                try
                    [U,W]=eig(X+Y2/mu2,'nobalance');
                catch
                    opts.tol=1e-4;
                    [U,W]=eigs(X+Y2/mu2,size(X,1)-1,'lm',opts);
                end
            end
            w=diag(W);
            ind=find(w>0);
            Z=U(:,ind)*diag(w(ind))*U(:,ind)';
        otherwise
            error(['lowRankPrior ' methodLowRankPrior ' not recognized'])
    end
            
    % update X
    X0 = X;
    X1 = (mu1*(M-E+Y1/mu1)+mu2*(Z-Y2/mu2))/(mu1+mu2);
    X2 = Z-Y2/mu2;
    for i = 1:n-1
        for j = i+1:n
            if flag(i,j)
                X(3*i-2:3*i,3*j-2:3*j) = X1(3*i-2:3*i,3*j-2:3*j);
            else
                X(3*i-2:3*i,3*j-2:3*j) = X2(3*i-2:3*i,3*j-2:3*j);
            end
            if innerProj
                X(3*i-2:3*i,3*j-2:3*j) = projectR(X1(3*i-2:3*i,3*j-2:3*j));
            end
            X(3*j-2:3*j,3*i-2:3*i) = X(3*i-2:3*i,3*j-2:3*j)';
        end
    end
    
    PrimRes1 = norm(S.*(M-X-E),'fro')/norm(M,'fro');
    PrimRes2 = norm(X-Z)/norm(M,'fro');
    DualRes = norm(X-X0,'fro')/norm(X0,'fro');
    if flagShowStats
        fprintf('Iter %d: Res = %f, %f, %f; mu = %f, %f; \n',iter,PrimRes1,PrimRes2,DualRes,mu1,mu2);
    end
    if flagCollectStats
        output.primRes(iter,:)=[PrimRes1,PrimRes2];
        output.dualRes(iter)=DualRes;
        output.mu(iter,:)=[mu1 mu2];
    end
    
    if  PrimRes1 < tol && PrimRes2 < tol && DualRes < tol
        break
    else
        if PrimRes1 > 10*DualRes
            mu1 = 2*mu1;
        elseif DualRes > 10*PrimRes1
            mu1 = mu1/2;
        end
        if PrimRes2 > 10*DualRes
            mu2 = 2*mu2;
        elseif DualRes > 10*PrimRes2
            mu2 = mu2/2;
        end
    end
    
end

if iter>=maxIter
    warning('Maximum number of iterations reached')
end

if innerProj == false
    for i = 1:n-1
        for j = i+1:n
            X(3*i-2:3*i,3*j-2:3*j) = projectR(X(3*i-2:3*i,3*j-2:3*j));
            X(3*j-2:3*j,3*i-2:3*i) = X(3*i-2:3*i,3*j-2:3*j)';
        end
    end
end

if flagCollectStats
    output.primRes(iter+1:end)=[];
    output.dualRes(iter+1:end)=[];
    output.mu(iter+1:end)=[];
    output.iter=iter;
end

end

function y = groupthresholdL2(x,lam)
nrm = norm(x(:));
y = max(nrm-lam,0)*x/(nrm+eps);
end
