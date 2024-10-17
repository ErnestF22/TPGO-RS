%function draw3dpolygon(V,F,color)
%Draws a 3D polyhedron having 3D vertices V (3 x Nvert array). The faces are triangle
%specified by the verteces indeces contained in F (3 x Nfaces array).
%COLOR is an optional MATLAB color specifier.

%%AUTORIGHTS%%

function h=draw3dpolygon(V,F,color,varargin)
flagwireframe=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'wireframe'
            flagwireframe=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

if(nargin<3)
    color='white';
end
NF=size(F,2);

X=[];
Y=[];
Z=[];


if(flagwireframe)
    for(iF=1:NF)
        X(:,iF)=V(1,F(:,iF));
        Y(:,iF)=V(2,F(:,iF));
        Z(:,iF)=V(3,F(:,iF));
    end
    h=plot3(X,Y,Z,color);
else
    for(iF=1:NF)
        for(v=1:3)
            X(v,iF)=V(1,F(v,iF));
            Y(v,iF)=V(2,F(v,iF));
            Z(v,iF)=V(3,F(v,iF));
        end
    end
    h=patch(X,Y,Z,color);
end
        