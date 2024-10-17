function [t,x,t_node,handleOde]=bnFranchi_Evolve(t_node,varargin)
TFinal=1;
odeSolver=@ode15s;
optsSolver={};
flagShowOdeProgress=false;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'tfinal'
            ivarargin=ivarargin+1;
            TFinal=varargin{ivarargin};
        case 'odesolver'
            ivarargin=ivarargin+1;
            odeSolver=varargin{ivarargin};
        case 'optssolver'
            ivarargin=ivarargin+1;
            optsSolver=[optsSolver varargin{ivarargin}];
        case 'showodeprogress'
            flagShowOdeProgress=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end
t_node.Ti0=t_node.Ti;
x0=t_node.Ti;

[d,NNodes]=size(x0);

%solve the closed loop ODE
if flagShowOdeProgress
    optsSolver = [optsSolver {'OutputFcn',@odeplot}];
end
 
handleOde=@(t,x) controller(t,x,t_node);
[t,x]=odeSolver(handleOde,[0 TFinal],x0(:),odeset(optsSolver{:}));
x=reshape(x',d,NNodes,[]);
t_node.Ti=x(:,:,end);

function dx=controller(t,x,t_node)
yd=t_node.Yijtruth;
d=size(yd,1);
gd=t_node.Gammaitruth;
E=t_node.E;
k=1;
x=reshape(x,d,[]);
NNodes=size(x,2);
[y,g]=bnFranchi_computeBearingAndGammaFormation(x);
y12=y(:,1);
y12d=yd(:,1);
dx=zeros(size(x));
dx(:,2)=k*orthComplementProjector(y12)*y12d;
for iNode=3:NNodes
    idxy=2+2*iNode-5;
    y1i=y(:,idxy);
    y1id=yd(:,idxy);
    dx(:,iNode)=k*(gd(iNode-2)*y1id-g(iNode-2)*y1i);
end
dx=reshape(dx,[],1);


