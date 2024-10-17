function [R,output]=sfm_rawRotationsLowRankALM_SolverL2(M,flag,varargin)
% this program solves
% min 1/2||M-R||^2 + lam ||R||_*
% st: each block of R corresponding to a non-zero block in M is in SO(3)
lambda = 10; % the regularization weight
tol=1e-4;
maxIter=3000;
flagShowStats=false;
flagCollectStats=nargout>1;
RInit=M;

output=[];

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        case 'tol'
            ivarargin=ivarargin+1;
            tol=varargin{ivarargin};
        case 'maxiter'
            ivarargin=ivarargin+1;
            maxIter=varargin{ivarargin};
        case 'init'
            ivarargin=ivarargin+1;
            RInit=varargin{ivarargin};
        case 'innerproj'
            ivarargin=ivarargin+1;
            innerProj=varargin{ivarargin};
        case 'showstats'
            flagShowStats=true;
        case 'collectstats'
            flagCollectStats=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

mu = .25/mean(abs(M(:)));
n=size(flag,1);

R = RInit;
X = R;
Y = zeros(size(M)); % the dual variable

if flagCollectStats
    output.maxIter=maxIter;
    output.primRes=NaN(maxIter,1);
    output.dualRes=NaN(maxIter,1);
    output.mu=NaN(maxIter,1);
end

% main
for iter = 1:maxIter
    
    % update R
    R0 = R;
    Z1 = (M+mu*X-Y)/(1+mu);
    Z2 = X-Y/mu; 
    for i = 1:n-1
        for j = i+1:n
            if flag(i,j)
                R(3*i-2:3*i,3*j-2:3*j) = Z1(3*i-2:3*i,3*j-2:3*j);
            else
                R(3*i-2:3*i,3*j-2:3*j) = Z2(3*i-2:3*i,3*j-2:3*j);
            end
            if innerProj
                R(3*i-2:3*i,3*j-2:3*j) = projectR(R(3*i-2:3*i,3*j-2:3*j));
            end
            R(3*j-2:3*j,3*i-2:3*i) = R(3*i-2:3*i,3*j-2:3*j)';
        end
    end
    
    % update X
    [U,W,V] = svd(R+Y/mu);
    w = diag(W) - lambda/mu;
    ind = find(w>0);
    X = U(:,ind)*diag(w(ind))*V(:,ind)';
    
    % update Y
    Y = Y + mu*(R-X);
    
    % Convergent?
    % We use the convergence criterion and parameter tuning scheme (for mu)
    % in section 3.4.1 of ADMM paper (Boyd 2010)
    primRes = norm(R-X,'fro')/norm(R0,'fro');
    dualRes = mu*norm(R-R0,'fro')/norm(R0,'fro');
    if flagShowStats
        fprintf('Iter %d: PrimRes = %f, DualRes = %f, rho = %f\n',iter,primRes,dualRes,mu);
    end
    
    if flagCollectStats
        output.primRes(iter)=primRes;
        output.dualRes(iter)=dualRes;
        output.mu(iter)=mu;
    end
    
    if  primRes < tol && dualRes < tol
        break
    else
        if primRes>10*dualRes
            mu = 2*mu;
        elseif dualRes>10*primRes
            mu = mu/2;
        else
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
