%Evolve trajectory of agents according to bearing control
%function [t,x,t_node,handleOde]=bearingNetworkEvolve(t_node, varargin)
%Inputs
%   t_node  structure with formation info
%Optional Inputs
%   'TFinal',t              maximum simulated time (default 1)
%   't_cost'                structure with information about cost to use
%       for the gradient, containing the fields
%       funsBearings        one of the structure of functions obtained with
%           bearingCostFunctions.m for the bearing terms in the cost
%       funsRanges          one of the structure of functions obtained with
%           bearingCostFunctions.m for the range terms in the cost
%           (optional)
%       flagUseRanges       boolean to enable/disable use of range
%           information (optional)
%       flagUseFeatures     boolean to enable/disable use of static features
%           information (optional, has effect only for dynamic case)
%   'odesolver', odexxx     function handle to one of the ODE solvers
%           (default ode15s)
%   'optsControl',cell      cell array of options to pass to
%           bearingNetworkControlDirect
%   'showODEProgress'       add 'OutputFcn',@odeplot to the options for the
%           ODE solver
%   'opstSolver'            cell array of options to pass to odeset
%TODO: allow range measurements
function [t,x,t_node,handleOde]=bearingNetworkEvolve(t_node, varargin)

flagZeroControl=false;

TFinal=1;
funsBearings=bearingCostFunctions('angleSq');
funsRanges=bearingCostFunctions('squared');
flagUseRanges=false;
optsControl={};
flagCollisionAvoidance=false;
alpha=1;

%options for the dynamic case
flagDynamic=false;
flagUseFeatures=false;
lambda=0;
mass=1;
optsControlDynamic={};

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
        case 't_cost'
            ivarargin=ivarargin+1;
            t_cost=varargin{ivarargin};
            flagUseRanges=false;
            flagUseFeatures=false;
            flagDynamic=false;
            
            funsBearings=t_cost.funsBearings;
            if isfield(t_cost,'funsRanges')
                funsRanges=t_cost.funsRanges;
                flagUseRanges=true;
            end
            if isfield(t_cost,'flagUseRanges')
                flagUseRanges=t_cost.flagUseRanges;
            end
            if isfield(t_cost,'flagUseFeatures')
                flagUseFeatures=t_cost.flagUseFeatures;
            end
            if isfield(t_cost,'flagDynamic')
                flagDynamic=t_cost.flagDynamic;
            end
            if isfield(t_cost,'flagCollisionAvoidance')
                flagCollisionAvoidance=t_cost.flagCollisionAvoidance;
            end
            if isfield(t_cost,'alpha')
                alpha=t_cost.alpha;
            end
        case 'optscontrol'
            ivarargin=ivarargin+1;
            optsControl=[optsControl varargin{ivarargin}];
        case 'optscontroldynamic'
            ivarargin=ivarargin+1;
            optsControlDynamic=[optsControlDynamic varargin{ivarargin}];
        case 'lambda'
            ivarargin=ivarargin+1;
            lambda=varargin{ivarargin};
        case 'mass'
            ivarargin=ivarargin+1;
            mass=varargin{ivarargin};
        case 'odesolver'
            ivarargin=ivarargin+1;
            odeSolver=varargin{ivarargin};
        case 'optssolver'
            ivarargin=ivarargin+1;
            optsSolver=[optsSolver varargin{ivarargin}];
        case 'showodeprogress'
            flagShowOdeProgress=true;
        case 'showargs'
            cellfun(@disp, varargin)
        case 'zerocontrol'
            flagZeroControl=true;
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

d=size(t_node.Ti,1);

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;

if flagUseRanges
    ERanges=t_node.Er;
    ygRanges=t_node.Yrijtruth;
    nygRanges=t_node.nYrijtruth;
end

if flagUseFeatures
    EFeatures=t_node.Ef;
    xFeatures=t_node.Tif;
end

if flagCollisionAvoidance
    ECollisionAvoidance=t_node.Eca;
    rCollisionAvoidance=t_node.rca;
end

optsControl=[optsControl 'alpha',alpha];

%initial conditions
t_node.Ti0=t_node.Ti;
x0=t_node.Ti;
if flagDynamic
    if isfield(t_node,'dTi')
        x0=[x0;t_node.dTi];
    else
        x0=[x0;zeros(size(x0))];
    end
end

%solve the closed loop ODE
if flagShowOdeProgress
    optsSolver = [optsSolver {'OutputFcn',@odeplot}];
end

%helper subfunction which implements the closed loop system
function dz=ode(t,z)
    %measurements
    if ~flagDynamic
        x=reshape(z,d,[]);
    else
        z=reshape(z,2*d,[]);
        x=z(1:d,:);
        dx=z(d+(1:d),:);
    end

    if flagZeroControl
        u=zeros(size(dx));
    else
        yBearings=bearingNetworkComputeBearings(x,EBearings);

        %control
        inpsControl={EBearings,yBearings,ygBearings,funsBearings};
        if flagUseRanges
            [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
            inpsControl=[inpsControl 'ranges',ERanges,yRanges,ygRanges,nyRanges,nygRanges,funsRanges];
        end
        u=bearingNetworkControlDirect(inpsControl{:},optsControl{:},'t',t);
        if flagCollisionAvoidance
            [yCA,aCA]=bearingNetworkCollisionImages(x,ECollisionAvoidance,rCollisionAvoidance);
            u=bearingNetworkCollisionProject(u,yCA,aCA,ECollisionAvoidance);
        end
        if flagDynamic
            dyBearing=measurementsBearingsDynamic(x,dx,EBearings);
            inpsControlDynamic={EBearings,dyBearing};
            if flagUseRanges
                dqRanges=bearingNetworkComputeRangeResidualsDerivatives(dx,ERanges,ygRanges);
                inpsControlDynamic=[inpsControlDynamic 'ranges',ERanges,dqRanges,ygRanges];
            end
            if flagUseFeatures
                dyFeatures=measurementsFeaturesBearings(x,dx,xFeatures,EFeatures);
                inpsControlDynamic=[inpsControlDynamic 'features',EFeatures,dyFeatures];
            end
            u=u+bearingNetworkControlDynamic(inpsControlDynamic{:},optsControlDynamic{:});
        end
    end
    if ~flagDynamic
        dz=u(:);
    else
        dz=[dx;u/mass-lambda*dx];
        dz=dz(:);
    end
end

%note: ode is the helper subfunction defined above
[t,x]=odeSolver(@ode,[0 TFinal],x0(:),odeset(optsSolver{:}));

if ~flagDynamic
    dState=d;
else
    dState=2*d;
end
x=reshape(x',dState,testNetworkGetNumberOfNodes(t_node),[]);

t_node.Ti=x(1:d,:,end);
handleOde=@ode;
end

%helper function for computing dynamic measurements of bearings
function dyBearing=measurementsBearingsDynamic(x,dx,EBearings)
[yBearings,nyBearing]=bearingNetworkComputeBearings(x,EBearings);
dyBearing=bearingNetworkComputeBearingsDerivative(dx,EBearings,yBearings,nyBearing);
end

%helper function for computing dynamic measurements of feature bearings
function dyFeatures=measurementsFeaturesBearings(x,dx,xFeatures,EFeatures)
[yFeatures,nyFeatures]=bearingNetworkComputeFeatureBearings(x,xFeatures,EFeatures);
dyFeatures=bearingNetworkComputeFeatureBearingsDerivative(dx,EFeatures,yFeatures,nyFeatures);
end
