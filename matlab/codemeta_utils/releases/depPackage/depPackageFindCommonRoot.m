%function root=depPackageFindCommonRoot(nameList)
%Given a cell array of strings, find the common root among all of them
function root=depPackageFindCommonRoot(nameList)

%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

% Ver: 09-Feb-2012 19:36:21

% Roberto Tron (tron@bu.edu)

%convert to char
charList=char(nameList);
%use diff to find how many columns have all the same character
rootLength=find(max(diff(double(charList))~=0),1,'first')-1;

root=nameList{1}(1:rootLength);
