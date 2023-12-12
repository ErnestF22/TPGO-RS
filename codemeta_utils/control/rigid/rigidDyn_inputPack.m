function dx=rigidDyn_inputPack(dw,dv,varargin)
flagExtendedVector=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'state'
            ivarargin=ivarargin+1;
            R=varargin{ivarargin};
            ivarargin=ivarargin+1;
            w=varargin{ivarargin};
            ivarargin=ivarargin+1;
            v=varargin{ivarargin};
            flagExtendedVector=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if ~flagExtendedVector
    dx=[dw;dv];
else
    dx=[rotDyn_statePack(multiprod(R,hat3(w)),dw);
        realDyn_statePack(v,dv)];
end
