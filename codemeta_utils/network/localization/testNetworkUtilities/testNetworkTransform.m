%function t_node=testNetworkTransform(t_node,gGlobal,direction)
%Applies gGlobal on the left or right (depending on the value of direction)
%to each one of the rigid body transformations in t_node.(member)
function t_node=testNetworkTransform(t_node,gGlobal,varargin)
structType=testNetworkDetectStructType(t_node);
N=testNetworkGetNumberOfNodes(t_node);

member='gi';
direction='left';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'member'
            ivarargin=ivarargin+1;
            member=lower(varargin{ivarargin});
        case 'direction'
            ivarargin=ivarargin+1;
            direction=lower(varargin{ivarargin});
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

for iNode=1:N
    switch structType
        case 'single'
            switch direction
                case 'left'
                    t_node.gi(:,:,iNode)=gGlobal*t_node.gi(:,:,iNode);
                case 'right'
                    t_node.gi(:,:,iNode)=t_node.gi(:,:,iNode)*gGlobal;
            end
        case 'array'
            switch direction
                case 'left'
                    t_node(iNode).(member)=gGlobal*t_node(iNode).(member);
                case 'right'
                    t_node(iNode).(member)=t_node(iNode).(member)*gGlobal;
            end
    end
end
