%Plot points with different styles according to group membership
%function plotGroups(x,idxX)
function plotGroups(x,idxX,varargin)
labels=unique(idxX);
NLabels=length(labels);
markerOrder={'o','x','+','d'};
NMarkers=length(markerOrder);
NColors=max(ceil(NLabels/NMarkers),5);
colorOrder=rbg(NColors);
styleOrder=cell(1,NColors*NMarkers);
cnt=1;
for iMarker=1:NMarkers
    for iColor=1:NColors
        styleOrder{cnt}={markerOrder{iMarker},'Color',colorOrder(iColor,:)};
        cnt=cnt+1;
    end
end

flagHold=ishold();
for iLabel=1:NLabels
    plotPoints(x(:,idxX==labels(iLabel)),styleOrder{iLabel}{:},varargin{:});
    hold on
end
if ~flagHold
    hold off
end

