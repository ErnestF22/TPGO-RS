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


function quickshift_checkMembershipsClusters(membershipCluster,membershipPrior,clusterIndicators)
%check that clusters with non-zero indicators are the same as
%the unique clusters referenced by the membership
clusters=unique(membershipCluster);
if exist('clusterIndicators','var')
    fprintf('Indeces of referred clusters...')
    clusterDiff=setdiff(clusters,find(sum(clusterIndicators,2)));
    if isempty(clusterDiff)
        disp('OK')
    else
        disp('Inconsistent')
    end
end

%check that each cluster does not contain points from the same prior
fprintf('Check that clusters contain points from different priors...')
flagPass=true;
for kCluster=clusters
    flagPoints=membershipCluster==kCluster;
    pointsPrior=membershipPrior(flagPoints);
    flagPass=and(flagPass,length(unique(pointsPrior))==length(pointsPrior));
end

if flagPass
    disp('OK')
else
    disp('Conflicts')
end
