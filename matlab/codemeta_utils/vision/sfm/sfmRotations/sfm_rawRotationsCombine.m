%Combine relative rotations into absolute rotations
%function [Ri,output]=sfm_rawRotationsCombine(Rij,E,varargin)
function [Ri,output]=sfm_rawRotationsCombine(Rij,E,varargin)
flagInferenceOutliers=false;
flagLocalRefine=false;
optsInference={};
optsGlobal={};
optsLocal={};

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'inferenceoutliers'
            flagInferenceOutliers=true;
        case 'flaginferenceoutliers'
            ivarargin=ivarargin+1;
            flagInferenceOutliers=varargin{ivarargin};
        case 'optsglobal'
            ivarargin=ivarargin+1;
            optsGlobal=[optsGlobal varargin{ivarargin}];
        case 'optsinference'
            ivarargin=ivarargin+1;
            optsInference=[optsInference varargin{ivarargin}];
        case 'localrefine'
            flagLocalRefine=true;
        case 'flaglocalrefine'
            ivarargin=ivarargin+1;
            flagLocalRefine=varargin{ivarargin};
        case 'optslocal'
            ivarargin=ivarargin+1;
            optsLocal=[optsInference varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if flagInferenceOutliers
    %make sure that edges are not repeated
    [E,idxE]=edges2edges(E,'toOriented');
    Rij=Rij(:,:,idxE);
    
    flagE=sfm_rawRotationsEdgeLabelsInference(Rij,E,optsInference{:});
    Rij=Rij(:,:,~flagE);
    E=E(~flagE,:);
    
    output.inference.flagE=flagE;
    output.inference.Rij=Rij;
    output.inference.E=E;
end

[Ri,output.factorization]=sfm_rawRotationsCombineFactorization(Rij,E,optsGlobal{:});

if flagLocalRefine
    output.RiFactorization=Ri;
    Ri=sfm_rawRotationsRefine(Rij,E,Ri,optsLocal{:});
end
