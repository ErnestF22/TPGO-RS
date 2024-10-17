%function c=medoids_centerCost(x,varargin)
%Compute the medoids cost for a given center, without assignments.
%Inputs
%   x   [d x nbPoints]  Data points
%Optional inputs
%   'bias',biasMult     [d x nbBias]    Multiplicative bias term (summed
%   with negative sign). If a matrix is passed, the columns are summed.
%   'weights',weights   [1 x nbPoints]  Weights for each absolute value term
%   'xEval',xEval       [d x nbPointsE] Points at which to evaluate the cost (default, x)
function c=medoids_centerCost(x,varargin)
nbPoints=size(x,2);
xEval=x;
flagBias=false;
flagPriorCenters=false;
weights=ones(1,nbPoints);
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'bias'
            ivarargin=ivarargin+1;
            biasMult=varargin{ivarargin};
            flagBias=true;
        case 'weights'
            ivarargin=ivarargin+1;
            weights=varargin{ivarargin};
        case 'priorcenters'
            ivarargin=ivarargin+1;
            priorCenters=varargin{ivarargin};
            ivarargin=ivarargin+1;
            priorCentersWeights=varargin{ivarargin};
            flagPriorCenters=true;
        case 'xeval'
            ivarargin=ivarargin+1;
            xEval=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if flagPriorCenters
    weights=[weights priorCentersWeights];
    x=[x priorCenters];
end

%compute cost from absolute value terms
c=sum(weights'.*distMatrixManhattan(x,xEval));
%compute bias
if flagBias
    biasMul=sum(biasMult,2);
    c=c-biasMul'*xEval;
end
