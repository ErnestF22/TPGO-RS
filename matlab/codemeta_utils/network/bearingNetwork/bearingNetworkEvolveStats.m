%Compute statistics from the trajectories provided by bearingNetworkEvolveStats
%function output=bearingNetworkEvolveStats(t_node,t,x,varargin)
%Inputs
%   t_node  structure for the formation
%   t,x     time and trajectory given by bearingNetworkEvolve
%Optional inputs
%   'cost'              compute the cost at each iteration
%   't_cost'                structure with information about cost to use
%       for the gradient, containing the fields
%       funsBearings        one of the structure of functions obtained with
%           bearingCostFunctions.m for the bearing terms in the cost
%       funsRanges          one of the structure of functions obtained with
%           bearingCostFunctions.m for the range terms in the cost
%           (optional)
%       flagUseRanges       boolean to enable/disable use of range
%           information (optional)
function output=bearingNetworkEvolveStats(t_node,t,x,varargin)
funsBearings=bearingCostFunctions('angleSq');
funsRanges=bearingCostFunctions('squared');
flagUseRanges=false;
alpha=1;

%options for the dynamic case
flagDynamic=false;

flagGetCost=false;
flagGetAngles=false;
flagGetConfigDistance=false;
flagGetResiduals=false;
flagGetRangeDifference=false;
flagGetRanges=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'cost'
            flagGetCost=true;
        case 't_cost'
            ivarargin=ivarargin+1;
            t_cost=varargin{ivarargin};
            flagUseRanges=false;
            flagDynamic=false;
            
            funsBearings=t_cost.funsBearings;
            if isfield(t_cost,'funsRanges')
                funsRanges=t_cost.funsRanges;
                flagUseRanges=true;
            end
            if isfield(t_cost,'flagUseRanges')
                flagUseRanges=t_cost.flagUseRanges;
            end
            if isfield(t_cost,'flagDynamic')
                flagDynamic=t_cost.flagDynamic;
            end
            if isfield(t_cost,'alpha')
                alpha=t_cost.alpha;
            end
        case 'angles'
            flagGetAngles=true;
        case 'bearingranges'
            flagGetRanges=true;
        case 'configurationdistance'
            flagGetConfigDistance=true;
        case 'residuals'
            flagGetResiduals=true;
        case 'rangedifference'
            flagGetRangeDifference=true;
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

EBearings=t_node.E;
ygBearings=t_node.Yijtruth;
nygBearings=t_node.nYijtruth;

if flagUseRanges
    ERanges=t_node.Er;
    ygRanges=t_node.Yrijtruth;
    nygRanges=t_node.nYrijtruth;
    NEdgesRanges=size(ERanges,1);
end

output.TFinal=max(t);
output.m=squeeze(sum(x,2))/t_node.NNodes;

NEdgesBearings=size(EBearings,1);
NNodes=t_node.NNodes;
Nit=length(t);
d=size(t_node.Ti,1);

if flagGetCost
    phi=zeros(Nit,1);
end
if flagGetAngles
    aEval=@(x) angles(x,EBearings,ygBearings);
    a=zeros(Nit,NEdgesBearings);
end
if flagGetResiduals
    c=zeros(Nit,NEdgesBearings);
    if flagUseRanges
        q=zeros(Nit,NEdgesRanges);
    end
end    
if flagGetConfigDistance
    distEval=@(x) bearingNetworkConfigurationDistance(t_node.Titruth,x,'flagUseRanges',flagUseRanges);
    dist=zeros(Nit,1);
end
if flagGetResiduals
    rd=zeros(Nit,NEdgesBearings);
end    
if flagGetRanges
    r=zeros(Nit,NEdgesBearings);
end    
%     w=getTextWaitBar(Nit);
%     w(0)

for it=1:Nit
    xIt=x(1:d,:,it);
    if flagDynamic
        dxIt=x(d+(1:d),:,it);
    end
        
    if flagGetCost
        phi(it)=cost(xIt);
        if flagDynamic
            phi(it)=phi(it)+0.5*sum(dxIt(:).^2);
        end            
        output.phi=phi;
    end
    if flagGetAngles
        a(it,:)=aEval(xIt);
        output.a=a;
    end
    if flagGetConfigDistance
        dist(it)=distEval(xIt);
        output.dist=dist;
    end
    if flagGetResiduals
        if ~flagUseRanges
            c(it,:)=residuals(x(:,:,it));
        else
            [c(it,:),q(it,:)]=residuals(x(:,:,it));
            output.q=q;
        end
        output.c=c;
    end
    if flagGetRangeDifference
        rd(it,:)=rangeDiff(xIt,EBearings,nygBearings);
        output.rd=rd;
    end
    if flagGetRanges
        r(it,:)=ranges(xIt);
        output.r=r;
    end
end

function c=cost(x)
    [yBearings,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
    if ~flagUseRanges
        c=bearingNetworkCost(EBearings,yBearings,ygBearings,nyBearings,funsBearings);
    else
        [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
        
        c=bearingNetworkCostCombined(...
            EBearings,ERanges,...
            yBearings,yRanges,ygBearings,ygRanges,...
            nyBearings,nyRanges,nygRanges,...
            funsBearings,funsRanges,alpha);
    end
    
end

function [c,q]=residuals(x)
    yBearings=bearingNetworkComputeBearings(x,EBearings);
    c=bearingNetworkComputeBearingsCosines(yBearings,ygBearings);

    if flagUseRanges
        [yRanges,nyRanges]=bearingNetworkComputeBearings(x,ERanges);
        q=bearingNetworkComputeRangeResiduals(yRanges,ygRanges,nyRanges,nygRanges);
    end
end

function r=ranges(x)
[~,r]=bearingNetworkComputeBearings(x,EBearings);
end

end

function a=angles(x,EBearings,ygBearings)
yBearings=bearingNetworkComputeBearings(x,EBearings);
a=sphere_dist(yBearings,ygBearings,'vector');
end

function rd=rangeDiff(x,EBearings,nygBearings)
[~,nyBearings]=bearingNetworkComputeBearings(x,EBearings);
rd=nyBearings-nygBearings;
end