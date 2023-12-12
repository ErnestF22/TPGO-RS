%Solve the low-rank least squares problem with ADMM and SVD projection
%Solves the following problem: find a matrix X by 
%min ||A*vec(X)-b||^2 
%subject to 
%   C*vec(X)>D
%   rank(X)=k
function [X,U,output]=fixedRankLSQP(A,b,C,d,sz,k,varargin)
flagDebugOutput=false;
if nargout>2
    flagDebugOutput=true;
end

flagProvidedReference=false;
flagProvidedInitial=false;
methodLS='overdetermined';

%parameters
rho=1;
maxIt=100;
tolInequalities=1e-6;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'referencesolution'
            ivarargin=ivarargin+1;
            XReference=varargin{ivarargin};
            flagProvidedReference=true;
        case 'initialsolution'
            ivarargin=ivarargin+1;
            XInitial=varargin{ivarargin};
            flagProvidedInitial=~isempty(XInitial);
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

nbInequalities=size(C,1);

%if a problem does not include inequalities, all parts of the code that
%regard these are effectively disabled by this flag
flagHasInequalities=nbInequalities>0;

%elements of the quadratic part from the LS problem
I=eye(size(A,2));
switch methodLS
    case 'compact'
        P=A'*A;
        q=A'*b;
    case 'overdetermined'
        AExtended=cat(1,A,repmat(rho*I,nbInequalities+1,1));
end

if flagHasInequalities
    %normalize inequalities
    [C,d]=linearInequalitiesNormalize(C,d);
end

%initialization
%X is the final solution
if ~flagProvidedInitial
    X=reshape(A\b,sz(1),sz(2));
else
    X=XInitial;
end
%U.FixedRank contains the Lagrange multipliers for the fixed rank constraint
U.FixedRank=zeros(sz);
%Z contains the projection of X onto the fixed rank set
Z=fixedRankProject(X,k);

if flagHasInequalities
    %UInequality contains the Lagrange multipliers for each inequality
    %constraints
    U.Inequalities=zeros(prod(sz),nbInequalities);
    %Y contains the projection of X on each one of the inequalities
    Y=linearInequalitiesProject(vec(X),C,d,'normalized');
end

%initialization of debug info
if flagDebugOutput
    output.linearResiduals=zeros(1,maxIt);
    output.fixedRankResiduals=zeros(1,maxIt);
    output.solutionSvals=zeros(min(sz),maxIt);
    if flagProvidedReference
        output.solutionReferenceResiduals=zeros(1,maxIt);
        output.projectedSolutionReferenceResiduals=zeros(1,maxIt);
    end
    if flagHasInequalities
        output.inequalitiesResiduals=zeros(nbInequalities,maxIt);
        output.inequalitiesNbActive=zeros(1,maxIt);
    end
end

%ADMM iterations
for it=1:maxIt
    vFixedRank=vec(Z-U.FixedRank);
    if flagHasInequalities
        flagInequalitiesActive=C*vec(X)-d<tolInequalities;
        vInequalities=Y-U.Inequalities;
        vInequalities(:,~flagInequalitiesActive)=Y(:,~flagInequalitiesActive);
    end
    switch methodLS
        case 'compact'
            if flagHasInequalities
                X=(P+(nbInequalities+1)*rho*I)\(rho*(vFixedRank+sum(vInequalities,2))+q);
            else
                X=(P+rho*I)\(rho*vFixedRank+q);
            end
        case 'overdetermined'
            bExtended=[b;vFixedRank];
            if flagHasInequalities
                bExtended=[bExtended; vec(vInequalities)];
            end
            X=AExtended\bExtended;
        otherwise
            error('methodLS not recognized')
    end
    X=reshape(X,sz(1),sz(2));
    Z=fixedRankProject(X+U.FixedRank,k);
    U.FixedRank=U.FixedRank+X-Z;
    if flagHasInequalities
        Y=linearInequalitiesProject(vec(X),C,d,'normalized');
        U.Inequalities=U.Inequalities+vec(X)-Y;
    end    
    %update debug info
    if flagDebugOutput
        output.linearResiduals(it)=norm(A*vec(X)-b,2)^2;
        output.fixedRankResiduals(it)=norm(X-Z,'fro');
        output.solutionSvals(:,it)=svd(X);
        if flagProvidedReference
            output.solutionReferenceResiduals(it)=norm(X-XReference,'fro');
            output.projectedSolutionReferenceResiduals(it)=norm(Z-XReference,'fro');
        end
        if flagHasInequalities
            output.inequalitiesResiduals(:,it)=min(C*vec(X)-d,0);
            output.inequalitiesNbActive(it)=sum(flagInequalitiesActive);
        end
    end
end

%project a matrix onto the set of rank(X)=k matrices
function Z=fixedRankProject(X,k)
%sanitize requested rank w.r.t. matrix dimensions
k=min(k,min(size(X)));
%compute projection using SVD
[U,S,V]=svd(X);
Z=U(:,1:k)*S(1:k,1:k)*V(:,1:k)';

