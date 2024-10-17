%Compute bearing control field for a single pair of nodes
%function dx=bearingNetworkFieldPair(y,yg,optsField)
function dx=bearingNetworkFieldPair(y,yg,varargin)
dx=[
    bearingField(y,yg,varargin{:})...
    bearingField(-y,-yg,varargin{:})...
   ];
