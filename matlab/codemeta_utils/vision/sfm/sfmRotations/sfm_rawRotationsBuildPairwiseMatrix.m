%Builds the matrices with relative rotations used for averargin
%function [GHat,D]=sfm_rawRotationsBuildPairwiseMatrix(Rij,E,varargin)
%Inputs
%   Rij     [3 x 3 x NEdges] vector of relative rotations
%   E       [NEdges x 2] list of edges
%Optional arguments
%   'NRotations',NRotations     Number of rotations generating Rij
%                               (default: max(E(:)))
%Outputs
%   GHat    [3*NRotations x 3*NRotation] matrix containing the Rij
%   D       [3*NRotations x 3*NRotation] diagonal matrix for the normalization
%See "Global Motion Estimation from Point Matches" for details.
function [GHat,D]=sfm_rawRotationsBuildPairwiseMatrix(Rij,E,varargin)
method='adjacency';
flagCountNormalize=true;
%Get the number of rotations from maximum number in edges
%Can be changed by passing corresponding argument
NRotations=max(E(:));

NEdges=size(E,1);

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nrotations'
            ivarargin=ivarargin+1;
            NRotations=varargin{ivarargin};
        case 'flagcountnormalize'
            ivarargin=ivarargin+1;
            flagCountNormalize=varargin{ivarargin};
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch method
    case 'adjacency'
        GHat = eye(3*NRotations);
    case 'laplacian'
        GHat = zeros(3*NRotations);
    otherwise
        error('method for building pairwise matrix not recognized')
end

count = zeros(NRotations);

idxRotation=reshape(1:3*NRotations,[3 NRotations]);
I=eye(3);

for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    iIdxRotation=idxRotation(:,iNode);
    jIdxRotation=idxRotation(:,jNode);
    
    switch method
        case 'adjacency'
            GHat(iIdxRotation,jIdxRotation)=GHat(iIdxRotation,jIdxRotation)+Rij(:,:,iEdge);
            GHat(jIdxRotation,iIdxRotation)=GHat(jIdxRotation,iIdxRotation)+Rij(:,:,iEdge)';
        case 'laplacian'
            GHat(iIdxRotation,jIdxRotation)=GHat(iIdxRotation,jIdxRotation)-Rij(:,:,iEdge);
            GHat(jIdxRotation,iIdxRotation)=GHat(jIdxRotation,iIdxRotation)-Rij(:,:,iEdge)';
            GHat(iIdxRotation,iIdxRotation)=GHat(iIdxRotation,iIdxRotation)+I;
            GHat(jIdxRotation,jIdxRotation)=GHat(jIdxRotation,jIdxRotation)+I;
    end
    count(iNode,jNode)=count(iNode,jNode)+1;
    count(jNode,iNode)=count(jNode,iNode)+1;
end

if flagCountNormalize
    for iNode=1:NRotations
        for jNode=find(count(iNode,:));
            iIdxRotation=idxRotation(:,iNode);
            jIdxRotation=idxRotation(:,jNode);
            
            GHat(iIdxRotation,jIdxRotation)=GHat(iIdxRotation,jIdxRotation)/count(iNode,jNode);
        end
    end
    count=count>0;
end

D=diag(kron(sum(count,2),ones(3,1)));