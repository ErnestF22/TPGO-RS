function output=relativePositionNetworkEvolveStats(t_node,t,x,varargin)
flagGetCost=false;
flagGetGradient=false;
flagGetControl=false;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'cost'
            flagGetCost=true;
        case 'gradient'
            flagGetGradient=true;
        case 'control'
            flagGetControl=true;
            ivarargin=ivarargin+1;
            handleOde=varargin{ivarargin};
        otherwise
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

E=t_node.E;
TijTruth=t_node.TijTruth;

output.TFinal=max(t);

    function c=cost(x)
        Tij=relativePositionNetworkCompute(x,E);
        c=relativePositionNetworkCost(E,Tij,TijTruth);
    end
if flagGetCost
    output.phi=funEvalVec(@cost,x);
end

    function g=gradient(x)
        Tij=relativePositionNetworkCompute(x,E);
        g=relativePositionNetworkCostGradient(E,Tij,TijTruth);
    end
if flagGetGradient
    output.gradPhi=funEvalVec(@gradient,x);
end

    function u=evalControl(t,x)
        u=handleOde(t,x(:));
        u=reshape(u,size(x));
    end
if flagGetControl
    Nit=size(x,3);
    control=zeros(size(x));
    for it=1:Nit
        control(:,:,it)=evalControl(t(it),x(:,:,it));
    end
    output.control=control;
end



end
