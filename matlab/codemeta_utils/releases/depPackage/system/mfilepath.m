%Return path of current mfile (usage similar to mfilename('fullpath'))
function p=mfilepath()

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

% Ver: 04-Jul-2018 12:49:05

% Roberto Tron (tron@bu.edu)

%get calling stack
s=dbstack('-completenames');
%s(1) is the current function
if length(s)==1
    %called from the command line
    p=pwd();
else
    %return only path of the calling file
    p=fileparts(s(2).file);
end


