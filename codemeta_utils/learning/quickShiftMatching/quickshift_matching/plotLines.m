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


%Plot lines given start and end points
%function plotLines(xStart,xEnd,varargin)
%This is essentially a smart version of plot which does not require
%separating and assembling the various coordinates of the points
function plotLines(xStart,xEnd,varargin)

if isempty(varargin) || ~isStyleString(varargin{1})
    %a style has not been provided, inject default marker and size
    varargin=[{'-'} varargin{:}];
elseif styleContainsLine(varargin{1})
    %a style has been provided, but it does not contain a line style
    %so add the default one
    varargin=[{['-' varargin{1}]} varargin{2:end}];
end

sz=size(xStart);
d=sz(1);
xData=zeros(2,prod(sz(2:end)),d);
for id=1:d;
    xData(:,:,id)=[squeeze(xStart(id,:,:));squeeze(xEnd(id,:,:))];
end
switch d
    case 2
        plot(xData(:,:,1),xData(:,:,2),varargin{:});
    case 3
        plot3(xData(:,:,1),xData(:,:,2),xData(:,:,3),varargin{:});
    otherwise
        error('First dimension of the data must be two or three')
end
