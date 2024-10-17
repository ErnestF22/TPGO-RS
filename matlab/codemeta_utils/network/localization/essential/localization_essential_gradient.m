function [t_node,output]=localization_essential_gradient(t_node,varargin)

optsLieMinimize={'disablePhase12'};
maxIt=10000;

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'displayit'
            optsLieMinimize=[optsLieMinimize 'displayIt'];
        case 'progressbar'
            optsLieMinimize=[optsLieMinimize 'progressBar'];
        case 'maxit'
            ivarargin=ivarargin+1;
            maxIt=varargin{ivarargin};
        case 'optslieminimize'
            ivarargin=ivarargin+1;
            optsLieMinimize=[optsLieMinimize varargin{ivarargin}];
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

E=testNetworkGetEdges(t_node);

[Ri0,Ti0]=testNetworkGetRotTransl(t_node);
GiInit=RT2G(Ri0,Ti0);

[~,Qij]=testNetworkGetRelativeEssential(t_node);

fEssential=@(G) essentialCostNetwork(G2R(G),G2T(G),Qij,E);

[GiEst,outputLieMinimize]=lie_minimizeGradNewton(rot3r3_funs(),fEssential,GiInit,...
    'epsilon',1,'MaxIt',maxIt,optsLieMinimize{:}); %,'showCost'

t_node.gi=GiEst;
output.outputLieMinimize=outputLieMinimize;

