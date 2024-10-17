%function draw3dcameraFromAxesAndCenter(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%   R   axes of the camera
%   T   center of the camera
%(this is the camera to world transformation, a.k.a. pose interpretation)
%Optional argument
%   'FlipY'         Flip the y axis when drawing the triangle
%   'Color1',c      Color for the sides of the pyramid
%   'Color2',c      Color for the base of the pyramid
%   'Scale',s       Scale for the pyramid (if 'Scale' is not specified, s=1)
%   'flagAxes',f    Flag to specify if axes should be drawn or not
%
%See also draw2dcameraFromPose

%%AUTORIGHTS%%

function draw2dcameraFromAxesAndCenter(R,T,varargin)
if(nargin<1)
    R=eye(2);
end
if(nargin<2)
    T=[0;0];
end

flagpostT=false;    %multiply T by R (corresponds to first rotate and then translate the camera)
flagFlipY=false;
flagAxes=true;
flagAxesLabels=true;
frustrumColor1=[15934	35723	14392]/65535/0.6;
scale=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'posttranslation'
            flagpostT=true;
        case 'color1'
            ivarargin=ivarargin+1;
           color1=varargin{ivarargin};
        case 'color2'
            ivarargin=ivarargin+1;
           color2=varargin{ivarargin};
        case 'flipy'
            flagFlipY=true;
        case 'scale'
            ivarargin=ivarargin+1;
            scale=varargin{ivarargin};
        case 'flagaxes'
            ivarargin=ivarargin+1;
            flagAxes=varargin{ivarargin};
        case 'flagaxeslabel'
            ivarargin=ivarargin+1;
            flagAxesLabels=varargin{ivarargin};
        case 'nolabels'
            flagAxesLabels=false;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

V=[0 1 -1;0 1 1];
F=[1 2 3; 2 3 1];

V=scale*V;

if flagFlipY
    V(2,:)=-V(2,:);
end

if(flagpostT)
    T=R*T;
end
V=R*V+T*ones(1,size(V,2));

flaghold=ishold;

h=draw2dpolygon(V,F);
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',frustrumColor1)
hold on
if flagAxes
    if flagAxesLabels
        plot2dframe(T,R)
    else
        plot2dframe(T,R,[],'nolabels')
    end
end

if(~flaghold)
    hold off
end
