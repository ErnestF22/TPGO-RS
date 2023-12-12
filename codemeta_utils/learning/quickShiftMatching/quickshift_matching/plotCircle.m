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


function varargout = plotCircle(x,y,r,varargin)
N=size(x,1);
if N>1
    h=repmat(struct(),N,1);
    flagHold=ishold;
    for iN=1:N
        plotCircle(x(iN,:),y(iN,:),r(iN),varargin{:})
        hold on
    end
    if ~flagHold
        hold off
    end
else
    d = r*2;
    px = x-r;
    py = y-r;
    h = rectangle('Position',[px py d d],'Curvature',[1,1],varargin{:});
end
if nargout>0
    varargout{1}=h;
end
