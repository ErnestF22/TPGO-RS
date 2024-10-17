function t_node=bearingNetworkAddMeasurements(t_node,varargin)
memberNode='Ti';
memberBearing='Yij';
memberEdges='E';

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
            case 'membernode'
                ivarargin=ivarargin+1;
                memberNode=varargin{ivarargin};
            case 'memberbearing'
                ivarargin=ivarargin+1;
                memberBearing=varargin{ivarargin};
            case 'memberedges'
                ivarargin=ivarargin+1;
                memberEdges=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end
E=t_node.(memberEdges);
Ti=t_node.(memberNode);

[Yij,nYij]=bearingNetworkComputeBearings(Ti,E);

t_node.(memberBearing)=Yij;
t_node.(['n' memberBearing])=nYij;
