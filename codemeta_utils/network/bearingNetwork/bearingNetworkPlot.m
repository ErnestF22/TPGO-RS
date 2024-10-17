%Plot a network together with its measurements
%function bearingNetworkPlot(t_node,varargin)
%Plot nodes using fields Titruth and (if available) Ti, together with edges
%between them indicating bearing measurements and (if available) range
%measurements.
%Optional arguments
%   'flagPlotRanges',flag       control plotting of range edges
%   'flagPlotBearings',flag     control plotting of bearing edges
%   'flagPlotFovs',flag         control plotting of fields of view
%           (default, false)

function bearingNetworkPlot(t_node,varargin)
flagPlotBearings=false;
flagPlotRanges=isfield(t_node,'Er');
flagPlotFeatureNodes=false;
flagPlotFeatureEdges=false;
flagPlotFovs=false;
flagLeader=false;
flagPlotCurrent=true;
flagPlotDesired=true;
colorCurrent=[1 0 0];
optsPlotNodes={};

%optional parameters
ivarargin=1;
while(ivarargin<=length(varargin))
    switch(lower(varargin{ivarargin}))
        case 'flagplotranges'
            ivarargin=ivarargin+1;
            flagPlotRanges=varargin{ivarargin};
        case 'flagplotfeaturenodes'
            ivarargin=ivarargin+1;
            flagPlotFeatureNodes=varargin{ivarargin};
        case 'flagplotfeatureedges'
            ivarargin=ivarargin+1;
            flagPlotFeatureEdges=varargin{ivarargin};
        case 'flagplotcurrent'
            ivarargin=ivarargin+1;
            flagPlotCurrent=varargin{ivarargin};
        case 'colorcurrent'
            ivarargin=ivarargin+1;
            colorCurrent=varargin{ivarargin};
        case 'flagplotdesired'
            ivarargin=ivarargin+1;
            flagPlotDesired=varargin{ivarargin};
        case 'plotfeatures'
            flagPlotFeatureEdges=true;
            flagPlotFeatureNodes=true;
        case 'flagplotbearings'
            ivarargin=ivarargin+1;
            flagPlotBearings=varargin{ivarargin};
        case 'flagplotfovs'
            ivarargin=ivarargin+1;
            flagPlotFovs=varargin{ivarargin};
        case 'leader'
            flagLeader=true;
            ivarargin=ivarargin+1;
            idxLeader=varargin{ivarargin};
        case 'optsplotnodes'
            ivarargin=ivarargin+1;
            optsPlotNodes=[optsPlotNodes varargin{ivarargin}];
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

E=t_node.E;
Titruth=t_node.Titruth;
Yijtruth=t_node.Yijtruth;

hold on

if flagPlotDesired
    bearingNetworkPlot_nodes(Titruth,'bs');
    if flagLeader && ~isempty(idxLeader)
        bearingNetworkPlot_nodes(Titruth(:,idxLeader),'bs','MarkerSize',5,'MarkerFaceColor','b')
    end
    bearingNetworkPlot_edges(Titruth,E,'b--')
    if isfield(t_node,'Er') && flagPlotRanges
        bearingNetworkPlot_edges(Titruth,t_node.Er,'b','linewidth',2)
    end

    if flagPlotBearings
        bearingNetworkPlot_bearing(Titruth,E,Yijtruth,'b')
    end
    if isfield(t_node,'Ritruth') && flagPlotFovs
        bearingPlotFov(Titruth,t_node.Ritruth,t_node.fov,'b',t_node.y0);
    end
else
    hold on
end

if isfield(t_node,'Tif') 
    if flagPlotFeatureNodes
        bearingNetworkPlot_nodes(t_node.Tif,'k');
    end
    if flagPlotFeatureEdges
        bearingNetworkPlot_features_edges(Titruth,t_node.Tif,t_node.Ef,'b:');
    end
end

if isfield(t_node,'Ti') && flagPlotCurrent
    bearingNetworkPlot_nodes(t_node.Ti,'o','Color',colorCurrent,optsPlotNodes{:})
    if flagLeader
        bearingNetworkPlot_nodes(t_node.Ti(:,idxLeader),'o','Color',colorCurrent,'MarkerSize',5,'MarkerFaceColor','r')
    end
    bearingNetworkPlot_edges(t_node.Ti,E,'-','Color',colorCurrent)
    if isfield(t_node,'Er') && flagPlotRanges
        bearingNetworkPlot_edges(t_node.Ti,t_node.Er,'-','linewidth',2,'Color',colorCurrent)
    end
    if flagPlotBearings
        bearingNetworkPlot_bearing(t_node.Ti,E,t_node.Yijtruth,'b')
        if isfield(t_node,'Yij')
            bearingNetworkPlot_bearing(t_node.Ti,E,t_node.Yij,'r')
        end
    end
    if isfield(t_node,'Ri') && flagPlotFovs
        bearingPlotFov(t_node.Ti,t_node.Ri,t_node.fov,'r',t_node.y0);
    end
    if isfield(t_node,'Tif') && flagPlotFeatureEdges
        bearingNetworkPlot_features_edges(t_node.Ti,t_node.Tif,t_node.Ef,'r:');
    end
end

hold off
axis equal

function bearingNetworkPlot_features_edges(Ti,Tif,E,varargin)
TijStart=Ti(:,E(:,1));
TijEnd=Tif(:,E(:,2));
plotLines(TijStart,TijEnd,varargin{:})

function bearingNetworkPlot_bearing(Ti,E,Yij,style)
TijStart=Ti(:,E(:,1));
quiver(TijStart(1,:),TijStart(2,:),Yij(1,:),Yij(2,:),0.5,style)
