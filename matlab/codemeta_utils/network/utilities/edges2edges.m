%Passes from symmetric to oriented edges (i.e., removes redundant edges)
%function [E,idxE]=edges2edges(E,varargin)
%Outputs
%   E       new edges
%   idxE    indeces of new edges in the old set
function [E,idxE]=edges2edges(E,varargin)
method='tooriented';

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'method'
            ivarargin=ivarargin+1;
            method=lower(varargin{ivarargin});
        case 'tooriented'
            ivarargin=ivarargin+1;
            method='tooriented';
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

idxE=(1:size(E,1))';
switch method
    case 'tooriented'
        iEdge=1;
        while iEdge<size(E,1)
            flagRepeated=E(iEdge+1:end,2)==E(iEdge,1) & E(iEdge+1:end,1)==E(iEdge,2);
            flagRepeated=flagRepeated | (E(iEdge+1:end,1)==E(iEdge,1) & E(iEdge+1:end,2)==E(iEdge,2));
            flagRepeated=[false(iEdge,1); flagRepeated];
            E(flagRepeated,:)=[];
            idxE(flagRepeated)=[];
            iEdge=iEdge+1;
        end
    otherwise
        error('edges2edges method not recognized')
end
