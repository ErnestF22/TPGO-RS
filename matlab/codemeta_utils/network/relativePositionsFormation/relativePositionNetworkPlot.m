function relativePositionNetworkPlot(t_node,varargin)
optsEdgesTi={};
optsEdgesTiTruth={};
flagEdgesTi=true;
flagEdgesTiTruth=true;

%optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'optsedgesti'
            ivarargin=ivarargin+1;
            optsEdgesTi=[optsEdgesTi varargin{ivarargin}];
        case 'optsedgestitruth'
            ivarargin=ivarargin+1;
            optsEdgesTiTruth=[optsEdgesTiTruth varargin{ivarargin}];
        case 'flagedgestitruth'
            ivarargin=ivarargin+1;
            flagEdgesTiTruth=varargin{ivarargin};
        case 'flagedgesti'
            ivarargin=ivarargin+1;
            flagEdgesTi=varargin{ivarargin};
        otherwise
            disp(varargin{ivarargin})
            error('Argument not valid!')
    end
    ivarargin=ivarargin+1;
end

if isempty(optsEdgesTi)
    optsEdgesTi={'r'};
end
if isempty(optsEdgesTiTruth)
    optsEdgesTiTruth={'b'};
end

if flagEdgesTi
    plotPoints(t_node.Ti,'r')
    plotEdgeArrows(t_node.Ti(:,:,end),t_node.E,optsEdgesTi{:});
end
if flagEdgesTiTruth
    plotPoints(t_node.TiTruth,'b');
    plotEdgeArrows(t_node.TiTruth,t_node.E,optsEdgesTiTruth{:});
end
