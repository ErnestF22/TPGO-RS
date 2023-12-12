function dx=realDyn_inputPack(dv,varargin)
flagExtendedVector=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'state'
            ivarargin=ivarargin+1;
            v=varargin{ivarargin};
            flagExtendedVector=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
if ~flagExtendedVector
    dx=dv;
else
    dx=realDyn_statePack(v,dv);
end
