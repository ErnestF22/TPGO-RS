%Compute distance between two configurations
%function d=bearingNetworkConfigurationDistance(xRef,x,varargin)
%Optional inputs
%   'flagUseRanges',flag
%       If false, compute distance in projective space
%       If true, compute Euclidean distance
function d=bearingNetworkConfigurationDistance(xRef,x,varargin)
flagUseRanges=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flaguseranges'
            ivarargin=ivarargin+1;
            flagUseRanges=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

NNodes=size(xRef,2);
xRef=xRef-mean(xRef,2)*ones(1,NNodes);
x=x-mean(x,2)*ones(1,NNodes);

if flagUseRanges
    d=norm(xRef-x,'fro');
else
    d=sphere_dist(cnormalize(xRef(:)),cnormalize(x(:)));
end
