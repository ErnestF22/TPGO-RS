%Compute bearing measurements and scales given node locations and incidence matrix
%function [u,lambda]=bearingCluster_getBearingsScalesFromB(x,B)
function [u,lambda]=bearingCluster_getBearingsScalesFromB(x,B,varargin)
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

%adjust sign of incident matrix
B=-B;

%augmented incidence matrix
Bd=kron(B,eye(d));

%location difference vectors along edges
rx=reshape(Bd*x(:),d,[]);

%add noise if requested
if flagNoisy
    rx=rx+sigmaNoise*randn(size(rx));
end

%compute bearings and scales
[u,lambda]=cnormalize(rx);
lambda=lambda';
