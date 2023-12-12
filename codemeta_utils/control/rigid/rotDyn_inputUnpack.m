function dw=rotDyn_inputUnpack(dx,varargin)
flagExtendedVector=size(dx,1)>3;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'flagextendedvector'
            ivarargin=ivarargin+1;
            flagExtendedVector=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if ~flagExtendedVector
    dw=dx(1:3,:);
else
    [~,dw]=rotDyn_stateUnpack(dx);
end

