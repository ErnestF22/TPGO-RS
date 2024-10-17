%function mu=medoids_centerSearchSingle(x,mu,kCenter,varargin)
%Computes the cost for all possible assignments of the k-th center, and
%assign to minimum.
%Inputs
%Optional inputs: Same as medoids.m
function [mu,cMin]=medoids_centersAssignmentSearchSingle(x,mu,kCenter,varargin)
muMin=mu;
cMin=medoids_centersAssignmentCost(x,muMin,varargin{:});
muPoint=mu;
xTrial=x;

%optional parameters to get prior centers
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'bias'
            ivarargin=ivarargin+1;
            %skip
        case 'priorcenters'
            ivarargin=ivarargin+1;
            priorCenters=varargin{ivarargin};
            xTrial=[xTrial [priorCenters{:}]];
            ivarargin=ivarargin+1;
            %skip weights
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end
nbPoints=size(xTrial,2);

for iPoint=1:nbPoints
    muPoint(:,kCenter)=xTrial(:,iPoint);
    cPoint=medoids_centersAssignmentCost(x,muPoint,varargin{:});
    if cPoint<cMin
        muMin=muPoint;
        cMin=cPoint;
    end
end
mu=muMin;
