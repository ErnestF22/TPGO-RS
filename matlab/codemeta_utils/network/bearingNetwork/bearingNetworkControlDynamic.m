%Compute velocity part of double integrator control law
%function ddx=bearingNetworkControlDynamic(EBearings,dyBearings,...
%   ERanges,dqRanges,ygRanges,alpha)
%Inputs
%   EBearings       edge list for bearing vectors measurements
%   dyBearings      derivative of bearing vectors corresponding to EBearings
%Optional Inputs
%   'Ranges',ERanges,dqRanges
%                   Add control from range measurements
%       ERanges         edge list for bearing vectors and range measurements
%       dqRanges        derivative of range residuals corresponding to ERanges
%       ygRanges        desired values for bearing vectors corresponding to ERanges
%   'Features',EFeatures,dyFeatures
%                   Add control from bearing to static features
%       EFeatures       edge list between nodes (indeces on first column)
%                       and features (indeces on second column)
%       dyFeatures      derivative of bearing vectors corresponding to
%                       EFeatures
%   'Alpha',alpha   [1 x 2] vector with weights for the two types of terms
%
%Output
%   ddx  control accelerations
%   
function ddx=bearingNetworkControlDynamic(EBearings,dyBearings,varargin)
flagUseRanges=false;
flagUseFeatures=false;
optsFeatures={};
optsBearings={};
alpha=1;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'ranges'
            flagUseRanges=true;
            ivarargin=ivarargin+1;
            ERanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dqRanges=varargin{ivarargin};
            ivarargin=ivarargin+1;
            ygRanges=varargin{ivarargin};
        case 'features'
            flagUseFeatures=true;
            ivarargin=ivarargin+1;
            EFeatures=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dyFeatures=varargin{ivarargin};
        case 'alpha'
            ivarargin=ivarargin+1;
            alpha=varargin{ivarargin};
        case 'optsfeatures'
            ivarargin=ivarargin+1;
            optsFeatures=[optsFeatures varargin{ivarargin}{:}];
        case 'optsbearings'
            ivarargin=ivarargin+1;
            optsBearings=[optsBearings varargin{ivarargin}{:}];
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if (flagUseRanges || flagUseFeatures) && length(alpha)==1
    alpha=[1 alpha];
end

if flagUseFeatures && length(alpha)==2
    alpha=[alpha 1];
end

NNodes=max(EBearings(:));
ddx=-alpha(1)*bearingNetworkControlDynamicBearings(EBearings,dyBearings,'NNodes',NNodes,optsBearings{:});
if flagUseRanges
    ddx=ddx-alpha(2)*bearingNetworkControlDynamicRanges(ERanges,dqRanges,ygRanges,'NNodes',NNodes);
end
if flagUseFeatures
    uFeatures=bearingNetworkControlDynamicFeatures(EFeatures,dyFeatures,'NNodes',NNodes,optsFeatures{:});
    ddx=ddx-alpha(3)*uFeatures;
end

