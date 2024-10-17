%Compute bearing measurements and scales given node locations and edges
%function [u,lambda]=bearingCluster_getBearingsScalesFromE(x,E)
function [u,lambda]=bearingCluster_getBearingsScalesFromE(x,E,varargin)
flagNoisy=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'noisy'
            ivarargin=ivarargin+1;
            flagNoisy=true;
            sigmaNoise=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%data dimension
d=size(x,1);

%location difference vectors along edges
rx=x(:,E(:,2))-x(:,E(:,1));

if flagNoisy
    rx=rx+sigmaNoise*randn(size(rx));
end

%compute bearings and scales
[u,lambda]=cnormalize(rx);
lambda=lambda';
