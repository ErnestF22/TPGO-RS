%Compute velocity part of double integrator control law from ranges
%function uRanges=bearingNetworkControlDynamicRanges(E,dq,yg,varargin)
function uRanges=bearingNetworkControlDynamicRanges(E,dq,yg,varargin)
NNodes=max(E(:));
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'nnodes'
            ivarargin=ivarargin+1;
            NNodes=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
d=size(yg,1);
uRanges=zeros(d,NNodes);
for iNode=1:NNodes
    iEdgeIn=E(:,1)==iNode;
    iEdgeOut=E(:,2)==iNode;
    uRanges(:,iNode)=-sum(yg(:,iEdgeIn).*(ones(d,1)*dq(iEdgeIn)'),2);
    uRanges(:,iNode)=uRanges(:,iNode)+sum(yg(:,iEdgeOut).*(ones(d,1)*dq(iEdgeOut)'),2);
end
