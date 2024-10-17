function [t,x,t_node,handleOde]=relativePositionNetworkEvolve(t_node,varargin)
TFinal=1;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'tfinal'
            ivarargin=ivarargin+1;
            TFinal=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

E=testNetworkGetEdges(t_node);
TijTruth=t_node.TijTruth;
Ti0=t_node.Ti;
[d,NNodes]=size(Ti0);

%prepare local function for closed loop 
    function dx=closedLoop(~,x)
        Ti=reshape(x,d,[]);
        Tij=relativePositionNetworkCompute(Ti,E);
        dTi=relativePositionNetworkControl(E,Tij,TijTruth);
        dx=dTi(:);
    end

handleOde=@closedLoop;


x0=Ti0(:);
odeSolver=@ode15s;
odeSolverOpts=odeset('OutputFcn',@odeplot);
[t,x]=odeSolver(@closedLoop,[0 TFinal],x0,odeSolverOpts);
x=reshape(x',d, NNodes, []);
t_node.Ti=x(:,:,end);

end