%Get approximated average absolute rotations from relative measurements
%function Ri=sfm_rawRotationsFactorize(G,varargin)
%Intput
%   G   [3*NRotations x 3*NRotations] array containing relative rotation
%       measurements
%Optional Inputs
%   'methodNormalization', method   Method for normalizing G before
%           extracting singular vectors
%       'none'      use the trailing singular vectors
%       'degree'    normalize with the degree of the node (computed as
%                   G(idxNode,:)*G(idxNode,:)')
%Output
%   Ri  [3 x 3 x NRotations] array of approximate rotations obtained from
%       the singular vectors
function Ri=sfm_rawRotationsFactorize(G,varargin)
methodNormalization='none';
methodSingularVectors='top';
flagGivenDNorm=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'methodnormalization'
            ivarargin=ivarargin+1;
            methodNormalization=lower(varargin{ivarargin});
        case 'methodsingularvectors'
            ivarargin=ivarargin+1;
            methodSingularVectors=lower(varargin{ivarargin});
        case 'dnorm'
            flagGivenDNorm=true;
            ivarargin=ivarargin+1;
            DNorm=lower(varargin{ivarargin});
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch methodNormalization
    case 'none'
        GNorm=G;
    case 'degree'
        if ~flagGivenDNorm
            DNorm=computeDegreeMatrix(G);
        end
        GNorm=DNorm\G/DNorm;
    case {'entropy','normalized'}
        [GNorm,output]=matProjectDoubleStochastic(G,methodNormalization);
        DNorm=output.DSqrtInv;
    otherwise
        error('methodFactorization not recognized')
end

GNorm=full(GNorm);
[~,S,V] = svd(GNorm);
%RiVec = V(:,1:3)*diag(sqrt(diag(S(1:3,1:3))));
switch methodSingularVectors
    case 'top'
        RiVec=V(:,1:3);
    case 'bottom'
        RiVec=V(:,end-2:end);
end

switch methodNormalization
    case 'none'
        %nothing to do
    case 'degree'
        RiVec=DNorm*RiVec;
end

Ri=sfm_rawRotationsDevectorize(RiVec);

function DNorm=computeDegreeMatrix(G)
NRotations=size(G,1);
idxRotations=reshape(1:3*NRotations,3,NRotations);
DNorm=zeros(3*NRotations);
for iRotation=1:NRotations
    idx=idxRotations(:,iRotation);
    DNorm(idx,idx)=sqrt(G(idx,:)*G(idx,:)');
end
