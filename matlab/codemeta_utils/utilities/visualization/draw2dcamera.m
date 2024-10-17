%function draw2dcamera(R,T,varargin)
%Draw a triangle and a set of axes to represent a camera
%Inputs
%   R  rotation of the camera
%   T  translation of the camera
function draw2dcamera(R,T,varargin)
if(nargin<1)
    R=eye(2);
end
if(nargin<2)
    T=[0;0];
end

flagpostT=false;    %multiply T by R (corresponds to first rotate and then translate the camera)

frustrumColor1=[15934	35723	14392]/65535;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'postTranslation'
            flagpostT=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

V=[0 1 -1;0 1 1];
F=[1 2 3; 2 3 1];

if(flagpostT)
    T=R*T;
end
V=R*V+T*ones(1,size(V,2));

flaghold=ishold;

h=draw2dpolygon(V,F);
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',frustrumColor1)
hold on
plot2dframe(T,R)

if(~flaghold)
    hold off
end

xlabel('x')
ylabel('y')

axis equal



