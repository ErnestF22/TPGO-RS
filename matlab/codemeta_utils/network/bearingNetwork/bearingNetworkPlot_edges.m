function bearingNetworkPlot_edges(Ti,E,varargin)
TijStart=Ti(:,E(:,1));
TijEnd=Ti(:,E(:,2));
plotLines(TijStart,TijEnd,varargin{:})
