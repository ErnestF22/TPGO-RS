%Infer outlier edges by sampling cycle closure errors
%function flagE=sfm_rawRotationsEdgeLabelsInference(R,E,varargin)
%The sampling of the cycle errors is obtained from
%sfm_rawRotationsSampleCycleErrors.m. The inference is done with
%sfm_rawRotationsCycleErrors2EdgeLabels.m.
%Inputs
%   R   [3 x 3 x NEdges] array of rotations for each edge
%   E   [NEdges x 2] or [MEdges x 3] list of (possibly weighted) edges
%Output
%   
%Optional inputs
%   'NSamples',N            Number of cycles (including triangles) to sample for
%                           doing the inference (default: 6*NEdges)
%   'optsInference', opts   Options to pass to inference step
function flagE=sfm_rawRotationsEdgeLabelsInference(R,E,varargin)
NSamples=6*size(E,1);
optsInference={};

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nsamples'
            ivarargin=ivarargin+1;
            NSamples=varargin{ivarargin};
        case 'optsinference'
            ivarargin=ivarargin+1;
            optsInference=[optsInference varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

%sample cycles
[C,e]=sfm_rawRotationsSampleCycleErrors(R,E,'triangles','NSamples',NSamples);
NCycles=size(C,2);

%inference
flagE=sfm_rawRotationsCycleErrors2EdgeLabels(C,e,optsInference{:});
