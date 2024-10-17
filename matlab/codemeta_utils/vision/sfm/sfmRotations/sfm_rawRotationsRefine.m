%Refine absolute rotations by minimizing sum of Riemannian distances
%function R=sfm_rawRotationsRefine(RRel,E,R,varargin)
%Inputs
%   RRel
%   E
%   R
%Optional inputs
%   'method',method     Type of cost to minimize
%       'weiszfeld'     Unsquared distances
%       'reshaped'      Distances reshaped with cost 'tron'
function R=sfm_rawRotationsRefine(RRel,E,R,varargin)
method='reshaped';
optsRefinement={'maxIt',1000,'RInit',R,'tolUpdate',1e-8};
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        case 'optsrefinement'
            ivarargin=ivarargin+1;
            optsRefinement=[optsInference varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

switch method
    case 'weiszfeld'
        R=rotLocWeiszfeldConsensus(E,RRel,optsRefinement{:});
    case 'reshaped'
        funs=rot3_almostGlobal_functions('type','tron','b',5);
        epsilon=1/(2*edges2maxDegree(E)*funs.mumax);
        R=rotLocRiemannianConsensus(E,RRel,optsRefinement{:},...
            'funs',funs,'epsilon',epsilon);
    otherwise
        error('Method for refinement not recognized')
end
