%Unpack the Euclidean and rotation parts from a vectorized form
%function [R,x]=rotdrd_stateUnpack(xVec,varargin)
%Optional inputs
%
function [R,x]=rotdrd_stateUnpack(xVec,varargin)
dReal=2;
nRot=1;
dRot=2;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'real'
            ivarargin=ivarargin+1;
            dReal=varargin{ivarargin};
        case 'rot'
            ivarargin=ivarargin+1;
            nRot=varargin{ivarargin};
            ivarargin=ivarargin+1;
            dRot=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

