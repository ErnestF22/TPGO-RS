% This file is part of QuickMatch.
%
% Copyright (C) 2018 Roberto Tron <tron@bu.edu> (Boston University)
% For more information see <https://bitbucket.org/tronroberto>
%
% QuickMatch is free software: you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation, either version 3 of the License, or
% (at your option) any later version.
%
% QuickMatch is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE. See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with QuickMatch. If not, see <http://www.gnu.org/licenses/>.


%Plot points with different styles according to group membership
%function plotGroups(x,idxX)
function plotGroups(x,idxX)
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
    plotPoints(x(:,idxX==labels(iLabel)),styleOrder{iLabel}{:});
    hold on
end
if ~flagHold
    hold off
end
