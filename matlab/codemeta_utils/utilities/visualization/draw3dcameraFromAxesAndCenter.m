%function draw3dcameraFromAxesAndCenter(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%R  axes of the camera
%T  center of the camera
%Optional argument
%   'FlipZ'         Flip the z axis when drawing the pyramid
%   'Color1',c      Color for the sides of the pyramid
%   'Color2',c      Color for the base of the pyramid
%   'Scale',s       Scale for the pyramid (if 'Scale' is not specified, s=1)
%   'flagAxes',f    Flag to specify if axes should be drawn or not
%
%See also draw3dcameraFromPose

%%AUTORIGHTS%%

function draw3dcameraFromAxesAndCenter(R,T,varargin)
if(nargin<1)
    R=eye(3);
end
if(nargin<2)
    T=[0;0;0];
end

flagpostT=false;    %multiply T by R (corresponds to first rotate and then translate the camera)
flagFlipZ=false;
flagAxes=true;
flagAxesLabels=true;
color1=[15934	35723	14392]/65535/0.6;
color2=[514	13107	39321]/65535;
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
        case 'flipz'
            flagFlipZ=true;
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

x=scale*[ 0  1  1 -1 -1];
y=scale*[ 0  1 -1 -1  1];
z=scale*[ 0  1  1  1  1];

if flagFlipZ
    z=-z;
end

F=[  1 2 3
     1 3 4;...
     1 4 5;...
     1 5 2]';

V=[x;y;z];
if(flagpostT)
    T=R*T;
end
V=R*V+T*ones(1,size(V,2));

flaghold=ishold;

if flagAxes
    if flagAxesLabels
        plot3dframe(T,R)
    else
        plot3dframe(T,R,[],'nolabels')
    end
end
hold on
h=draw3dpolygon(V,F);
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',color1)
h=patch(V(1,2:end),V(2,2:end),V(3,2:end), 'b');
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',color2)

if(~flaghold)
    hold off
end
