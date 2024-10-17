function draw3dCameraLineBackprojection(R,T,l2d,varargin)
optsDraw3dPlane={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'optsdraw3dplane'
            ivarargin=ivarargin+1;
            optsDraw3dPlane=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

%pose interpretation
l3d=[R'*l2d;T'*l2d];
o=homogeneousProjectPoints(R'*([0;0;1]-T),l3d);
draw3dPlane(l3d,'center',o,optsDraw3dPlane{:})
