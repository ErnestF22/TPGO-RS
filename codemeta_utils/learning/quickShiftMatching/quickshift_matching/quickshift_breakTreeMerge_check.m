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


function flag=quickshift_breakTreeMerge_check(treeComponents,componentData,componentIndicator)
flag=true;
if any(cellfun(@isempty,{componentData{treeComponents}}))
    disp('Checking that treeComponents points to non-empty componentData')
    disp('Error: invalid references to empty componentData')
    flag=false;
end

if ~all(cellfun(@isempty,componentData)==cellfun(@isempty,componentIndicator))
    disp('Checking that componentData and componentIndicator have same support')
    disp('Error: different supports (the two structures are not consistent)')
    flag=false;
end
