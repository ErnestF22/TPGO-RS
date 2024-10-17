function vOrth=triang_tangentBasisOrth(x,varargin)
flagNormalize=true;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagnormalize'
            ivarargin=ivarargin+1;
            flagNormalize=varargin{ivarargin};
        otherwise
                error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

vOrth=([-x(4); x(3); x(2); -x(1)]);
if flagNormalize
    vOrth=cnormalize(vOrth);
end

