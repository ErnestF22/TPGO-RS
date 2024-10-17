%function [C,E]=adj2incmatrix(A,varargin)
%Computes the directed [NEdges x NNodes] incidence matrix of a graph given the adjacency
%matrix. Each row contains a "1" for the first endpoint and a "-1" for the
%second endpoint for each edge (the graph is assumed to be directed).
%
%Optional parameters
%   'undirected'    Assume undirected graph. Look only at the upper
%                   triangular part of A, and do not include signs.
%   'oriented'      Assume undirected graph, but include signs.
%   'directed'      Assume directed graph.
function [C,E]=adj2incmatrix(A,varargin)
methodDirection='directed';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case {'undirected','oriented','directed'}
            methodDirection=lower(varargin{ivarargin});
        otherwise
            error('MATLAB:ArgumentInvalid',['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

E=adj2edges(A,methodDirection);

C=edges2incmatrix(E,size(A,1),methodDirection);

