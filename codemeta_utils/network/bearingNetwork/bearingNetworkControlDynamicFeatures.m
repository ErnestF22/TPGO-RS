%Compute velocity part of double integrator control law from feature bearings
%function uBearings=bearingNetworkControlDynamicFeatures(E,dy,varargin)
function uBearings=bearingNetworkControlDynamicFeatures(E,dy,varargin)
threshold=Inf;
flagClamp=false;
flagRemove=false;
flagDebug=false;

NNodes=max(E(:));
%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
        case 'threshold'
            ivarargin=ivarargin+1;
            threshold=varargin{ivarargin};
        case 'clamp'
            flagClamp=true;
        case 'remove'
            flagRemove=true;
        case 'debug'
            flagDebug=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
d=size(dy,1);
uBearings=zeros(d,NNodes);

if flagClamp
    dy=cnormalizeClamp(dy,threshold);
end
if flagRemove
    ndy=sqrt(sum(dy.^2));
    flagRemove=ndy>threshold;
    dy(:,flagRemove)=0;
    if flagDebug && sum(flagRemove)>0
        disp('Features removed:')
        disp(find(flagRemove))
    end
end

for iNode=1:NNodes
    uBearings(:,iNode)=-sum(dy(:,E(:,1)==iNode),2);
end
