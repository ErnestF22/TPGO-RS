function [GEst,output]=sfm_rawRotationsLowRankALM(G,varargin)

% parameters:
% - criterion: 'l2' (default), 'l1' or 'l12', corresponding to squared
% loss, L1 and group sparsity, respectively
% - innerProj: if doing SO(3) projection in each iteration
% - cvxInit: if using convex solution as initialization

criterion='l2';
innerProj=true;
cvxInit=true;
optsOptimizer={};

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'criterion'
            ivarargin=ivarargin+1;
            criterion=lower(varargin{ivarargin});
        case 'innerproj'
            ivarargin=ivarargin+1;
            innerProj=varargin{ivarargin};
        case 'cvxinit'
            ivarargin=ivarargin+1;
            cvxInit=varargin{ivarargin};
        case 'optsoptimizer'
            ivarargin=ivarargin+1;
            optsOptimizer= varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

flag=computeOccupancy(G);
switch criterion
    case {'l1','l12'}
        optsOptimizer = [optsOptimizer 'groupSparse' strcmp(criterion,'l12')];
        if innerProj
            if cvxInit 
                [GEst,output.solverInit]=sfm_rawRotationsLowRankALM_SolverL1(G,flag,optsOptimizer{:},'tol',1e-3,'innerProj',false);
                [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL1(G,flag,optsOptimizer{:},'innerProj',true,'init',GEst);
            else
                [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL1(G,flag,optsOptimizer{:},'innerProj',true);
            end
        else
            [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL1(G,flag,optsOptimizer{:},'innerProj',false);
        end
    case 'l2'
        if innerProj
            if cvxInit 
                [GEst,output.solverInit]=sfm_rawRotationsLowRankALM_SolverL2(G,flag,optsOptimizer{:},'tol',1e-3,'innerProj',false);
                [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL2(G,flag,optsOptimizer{:},'innerProj',true,'init',GEst);
            else
                [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL2(G,flag,optsOptimizer{:},'innerProj',true);
            end
        else
            [GEst,output.solver]=sfm_rawRotationsLowRankALM_SolverL2(G,flag,optsOptimizer{:},'innerProj',false);
        end
    otherwise
        error('criterion should be l2, l1 or l12.');
end

function flag=computeOccupancy(G)
NRotations=size(G,1)/3;
idxRotations=reshape(1:3*NRotations,3,NRotations);
flag=zeros(NRotations);
for iRotation=1:NRotations
    for jRotation=1:NRotations
        idxI=idxRotations(:,iRotation);
        idxJ=idxRotations(:,jRotation);
        if norm(G(idxI,idxJ),'fro')>eps
            flag(iRotation,jRotation)=1;
        end
    end
end

