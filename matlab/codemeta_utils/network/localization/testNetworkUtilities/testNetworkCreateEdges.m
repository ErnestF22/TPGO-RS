%function [E]=testNetworkCreateEdges(A,varargin)
%Creates the [NEdges x 2] matrix containing the endpoints of the edges in A
%Self-loops are discarded
%
%Optional arguments
%   'undirected'    Add only one direction for each edge. E.g. [1 2] but
%                   not [2 1]
%
function [E]=testNetworkCreateEdges(A,varargin)
flagUndirected=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'undirected'
            flagUndirected=true;
        case 'directed'
            flagUndirected=false;
        otherwise
            error('MATLAB:ArgumentInvalid',['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

[I,J]=find(A~=0);
if ~flagUndirected
    idxValid=I~=J;
else
    idxValid=I<J;
end
I=I(idxValid);
J=J(idxValid);

E=[I J];
