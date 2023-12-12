function [problem,muClusters]=admmMed_dataset(varargin)
nbNodes=4;
nbDatapoints=10;
nbClusters=5;
dimDatapoints=2;
sigmaDatapoints=0.02;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nbnodes'
            ivarargin=ivarargin+1;
            nbNodes=varargin{ivarargin};
        case 'nbdatapoints'
            ivarargin=ivarargin+1;
            nbDatapoints=varargin{ivarargin};
        case 'nbclusters'
            ivarargin=ivarargin+1;
            nbClusters=varargin{ivarargin};
        case 'dimdatapoints'
            ivarargin=ivarargin+1;
            dimDatapoints=varargin{ivarargin};

        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

muClusters=rand(dimDatapoints,nbClusters);

X=cell(1,nbNodes);
for iNode=1:nbNodes
    idxClusters=randi(nbClusters,1,nbDatapoints);
    X{iNode}=sigmaDatapoints*randn(dimDatapoints,nbDatapoints)+muClusters(:,idxClusters);
end

problem.nbNodes=nbNodes;
problem.nbDatapoints=nbDatapoints*ones(1,nbClusters);
problem.nbClusters=nbClusters;
problem.dimDatapoints=dimDatapoints;
problem.X=X;
