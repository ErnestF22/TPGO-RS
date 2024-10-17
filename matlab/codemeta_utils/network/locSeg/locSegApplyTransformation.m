%Apply given tranformation to the nodes of a network
%function [x,R]=locSegApplyTransformation(x,R,type,parameters,nodeList)
%Inputs
%   type        string containing the type of transformation to apply. Can
%               be one of the following
%       'rotation'
%       'translation'
%       'scale'
%   parameters  cell array containing parameteres for the transformation
%       'rotation'
%           parameters{1}   [3 x 1] rotation angle/axis representation
%           paramteres{2}   [3 x 1] rotation center
%       'translation'
%           parameters{1}   [3 x 1] translation vector
%       'scale'
%           parameters{1}   [1 x 1] scaling factor
%           parameters{2}   [3 x 1] scaling center
%   nodeList    row array with list of nodes affected by the transformation
%               (default: all the nodes)
function [x,R]=locSegApplyTransformation(x,R,type,parameters,nodeList)
if ~exist('nodeList','var')
    NNodes=size(x,2);
    nodeList=1:NNodes;
end

switch lower(type)
    case 'rotation'
        v=parameters{1};
        c=parameters{2};
        RTransformation=rot(v);
        for iNode=nodeList
            x(:,iNode)=RTransformation*(x(:,iNode)-c)+c;
            R(:,:,iNode)=RTransformation*R(:,:,iNode);
        end
    case 'translation'
        T=parameters{1};
        for iNode=nodeList
            x(:,iNode)=x(:,iNode)+T;
        end
    case 'scale'
        a=parameters{1};
        c=parameters{2};
        for iNode=nodeList
            x(:,iNode)=a*(x(:,iNode)-c)+c;
        end
    otherwise
        error('Type of requested transformation not valid')
end
        