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


function [flag, componentDataMerged]=quickshift_breakTreeMerge_fPrior(cmp1,cmp2,edge1,edge2)
%if a component is empty, initialize with the corresponding edge value
if isempty(cmp1)
    cmp1=edge1;
end
if isempty(cmp2)
    cmp2=edge2;
end

%flag=isempty(intersect(cmp1,cmp2));
flag=~overlapInt(cmp1,cmp2);
if flag
    componentDataMerged=[cmp1 cmp2];
else
    componentDataMerged=NaN;
end

%Return true if the intersection between two sets is non-empty
%The two sets are assumed to contain only positive integers
%Essentially, check the sum between two sparse indicator vectors to look
%for overlaps.
function flag=overlapInt(cmp1,cmp2)
l=max([cmp1 cmp2]);
l1=length(cmp1);
l2=length(cmp2);
u1=ones(1,l1);
u2=ones(1,l2);
flag=any((sparse(u1,cmp1,u1,1,l)+sparse(u2,cmp2,u2,1,l))>1);
