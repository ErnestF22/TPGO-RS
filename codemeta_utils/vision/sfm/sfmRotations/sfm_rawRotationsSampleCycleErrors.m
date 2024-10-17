%Samples cycles and computes closure errors
%function [C,e]=sfm_rawRotationsSampleCycleErrors(R,E)
%Samples cycles in the network by using Minimum Spanning Trees (MSTs). The
%weights for the first MST can be specified. The weights for the following
%MSTs are given by the number of times each edge has appeared in previously
%sampled cycles. The function also computes the closure error for the
%rotations along each sample cycle.
%Inputs
%   R   [d x d x NEdges] array of pairwise rotations
%   E   [NEdges x 2] or [NEdges x 3] list of edges. If the third column is
%       present, it indicates the weight for the first MST
%Optional Inputs
%   'NSamples',N    Maximum number of cycles to sample
%Outputs
%   C   [NEdges x NSamples] matrix with cycle vectors of sampled cycles
%   e   [NSamples x 1] vector of closure errors (radiants)
function [C,e]=sfm_rawRotationsSampleCycleErrors(R,E,varargin)
NSamples=100;
flagTriangles=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nsamples'
            ivarargin=ivarargin+1;
            NSamples=lower(varargin{ivarargin});
        case 'triangles'
            flagTriangles=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

C=grCycleSamples(E,'NSamples',NSamples);
if flagTriangles
    ETriangles=edges2triangles(E);
    C=[triangles2indicatorVectors(E,ETriangles) C];
end

NSamples=size(C,2);
l=zeros(NSamples,1);
e=zeros(NSamples,1);
for iSample=1:NSamples
    c=C(:,iSample);
    [idxc,signc]=edgesSortIndex(E,c);
    RCombined=sfm_rawRotationsCycleCombine_IndexSign(R,idxc,signc);
    l(iSample)=sum(abs(c));
    e(iSample)=rot_dist(RCombined,eye(3));
end

