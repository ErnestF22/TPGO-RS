%Transform rotation and translation pairs into rigid transformation matrices
%function G=RT2G(R,T,varargin)
%Inputs
%   R   [3 x 3 x N] array of rotation matrices
%   T   [3 x N] array of translation vectors
%Output
%   G   [4 x 4 x N] array of rigid trasformation matrices
%Optional Inputs
%   'compact'   Output G is a [3 x 4 x N] array
function G=RT2G(R,T,varargin)
flagCompact=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'compact'
            flagCompact=true;
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

N=size(R,3);
d=size(R,1);
G=[R permute(T,[1 3 2])];
if ~flagCompact
    G=[G; zeros(1,d,N) ones(1,1,N)];
end

