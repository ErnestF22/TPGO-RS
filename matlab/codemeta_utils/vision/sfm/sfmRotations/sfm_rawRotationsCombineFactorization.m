%Get stacked approximated average absolute rotations from relative measurements
%function [Ri,output]=sfm_rawRotationsCombineFactorization(Rij,E,varargin)
%Returns vector of rotations.
%Inputs
%   Rij     [3 x 3 x NEdgeRotations] array with relative rotations
%   E       [NEdgeRotations x 2] list of edges
%Optional Inputs
%   'methodLowRank', method     Method to use to improve the measurement
%           matrix
%       'Spectral'              no method applied (use directly the SVD step)
%       'SDP'                   apply semi-definite constraint, i.e.,
%                               sfm_rawRotationsLowRankSDP()
%       'ALM'                   Aumented Lagrange Multipliers method
%   'methodLowRankOpts', opts   cell array of options to pass to the method
%                               above (for now, only applies to 'ALM')
%   'optsFactorization', opts   cell array of options passed to
%                               sfm_rawRotationsFactorize
%   'normalize'                 add {{'methodNormalization','degree'}} to
%                               options for factorization
%   'adjust'                    use sfm_rawRotationsAdjust
%   'nrotations', NRotations    Number of absolute rotations to estimate 
%                               (if not provided, take N=max(E(:))
%   'projectBeforeFactorization'    project on SO(3) each 3x3 block of the
%                                   matrix produced by the method above
%Output
%   Ri   [3 x 3 x NRotations] absolute rotations ('reference' interpretation)
%
%See also sfm_rawRotationsLowRankSDP
function [Ri,output]=sfm_rawRotationsCombineFactorization(Rij,E,varargin)
flagGivenGHatAndD=false;
flagAdjust=false;
flagProjectBeforeFactorization=false;

methodLowRank='none';
optsLowRank={};
optsFactorization={};
flagCollectOutput=nargout>1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'ghatandd'
            ivarargin=ivarargin+1;
            G=varargin{ivarargin};
            ivarargin=ivarargin+1;
            D=varargin{ivarargin};
            flagGivenGHatAndD=true;
        case 'methodlowrank'
            ivarargin=ivarargin+1;
            methodLowRank=lower(varargin{ivarargin});
        case {'sdp','sdp2','alm','laplacian'}
            methodLowRank=lower(varargin{ivarargin});
        case 'spectral'
            methodLowRank='none';
        case 'projectbeforefactorization'
            flagProjectBeforeFactorization=true;
        case 'optslowrank'
            ivarargin=ivarargin+1;
            optsLowRank=[optsLowRank varargin{ivarargin}]; %#ok<AGROW>
        case 'optsfactorization'
            ivarargin=ivarargin+1;
            optsFactorization=[optsFactorization varargin{ivarargin}]; %#ok<AGROW>
        case 'normalize'
            optsFactorization=[optsFactorization 'methodNormalization' 'degree']; %#ok<AGROW>
        case 'adjust'
            flagAdjust=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%For using the method 'laplacian', we need a different pairwise matrix, and
%to consider the bottom singular vectors instead of the top
switch methodLowRank
    case 'laplacian'
        optsBuildPairwiseMatrix={'method','laplacian'};
        optsFactorization=[optsFactorization 'methodSingularVectors' 'bottom'];
    otherwise
        optsBuildPairwiseMatrix={};
end

if ~flagGivenGHatAndD
    [G,D]=sfm_rawRotationsBuildPairwiseMatrix(Rij,E,optsBuildPairwiseMatrix{:});
    optsFactorization=[optsFactorization, 'DNorm', D];
    if flagCollectOutput
        output.G=G;
        output.D=D;
    end
end

switch methodLowRank
    case {'spectral','none','laplacian'}
        GEst=G;
    case 'sdp'
        GEst=sfm_rawRotationsLowRankSDP(G,optsLowRank{:});
    case 'alm'
        [GEst,output.alm]=sfm_rawRotationsLowRankALM(G,optsLowRank{:});
    otherwise
        error('methodLowRank not recognized')
end

if flagCollectOutput
    output.GEst=GEst;
end

if ~flagProjectBeforeFactorization
    GEstProjected=GEst;
else
    GEstProjected=rotBlockProj(GEst);
end

RiFactorization=sfm_rawRotationsFactorize(GEstProjected,optsFactorization{:});

if flagAdjust
    RiAdjusted=sfm_rawRotationsAdjust(RiFactorization);
else
    RiAdjusted=RiFactorization;
end

Ri=projectR(RiAdjusted,'positiveDeterminant');

if flagCollectOutput
    output.RiFactorization=RiFactorization;
    output.RiAdjusted=RiAdjusted;
end

function GProj=rotBlockProj(G)
NRotations=size(G,1)/3;
GProj=zeros(size(G));
idxRotation=reshape(1:3*NRotations,[3 NRotations]);
for iRotation=1:NRotations
    iidx=idxRotation(:,iRotation);
    GProj(iidx,iidx)=rot_proj(G(iidx,iidx));
    for jRotation=iRotation+1:NRotations
        jidx=idxRotation(:,jRotation);
        GProj(iidx,jidx)=rot_proj(G(iidx,jidx));
        GProj(jidx,iidx)=GProj(iidx,jidx)';
    end
end


