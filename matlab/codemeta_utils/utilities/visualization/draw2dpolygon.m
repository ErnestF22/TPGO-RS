%function draw2dpolygon(V,F,color)
%Draws a 2-D polygon having 2-D vertices V (2 x Nvert array). The faces are lines
%specified by the verteces indeces contained in F (2 x Nfaces array).
%Other inputs
%   color   an optional MATLAB color specifier.

%%AUTORIGHTS%%

function h=draw2dpolygon(V,F,color,varargin)
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
    for iF=1:NF
        X(:,iF)=V(1,F(:,iF));
        Y(:,iF)=V(2,F(:,iF));
    end
    h=plot(X,Y,color);
else
    for iF=1:NF
        for v=1:2
            X(v,iF)=V(1,F(v,iF));
            Y(v,iF)=V(2,F(v,iF));
        end
    end
    h=patch(X',Y',color);
end
