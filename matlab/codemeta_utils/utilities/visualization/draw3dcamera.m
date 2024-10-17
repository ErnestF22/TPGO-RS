%function draw3dcamera(R,T,varargin)
%Draw a pyramid and a set of axes to represent a camera
%R  rotation of the camera
%T  translation of the camera

function draw3dcamera(R,T,varargin)
if(nargin<1)
    R=eye(3);
end
if(nargin<2)
    T=[0;0;0];
end

flagpostT=false;    %multiply T by R (corresponds to first rotate and then translate the camera)

frustrumColor1=[15934	35723	14392]/65535;
frustrumColor2=[514	13107	39321]/65535;
%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(varargin{ivarargin})
        case 'postTranslation'
            flagpostT=true;
%         case 'x0'
%             ivarargin=ivarargin+1;
%             x0=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

x=[ 0  1  1 -1 -1];
y=[ 0  1 -1 -1  1];
z=[ 0  1  1  1  1];

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

plot3dframe(T,R)
hold on
h=draw3dpolygon(V,F);
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',frustrumColor1)
h=patch(V(1,2:end),V(2,2:end),V(3,2:end),'b');
set(h,'FaceAlpha',0.5)
set(h,'FaceColor',frustrumColor2)

if(~flaghold)
    hold off
end

xlabel('x')
ylabel('y')
zlabel('z')

axis equal
view(3)



