%An example that packages depPackage itself
function depPackage_test

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

% Ver: 09-Jul-2018 12:53:17

% Roberto Tron (tron@bu.edu)

%name of the package
opts.packageName='depPackage';
%define directory remappings
opts.pairs={...
    'depPackage/','./'...
    };
%staging directory relative to where the current script is
opts.baseDir=fullfile(mfilepath(),'..','..','releases','depPackage');
opts.autorightsOpts={'dateVersion','authorName','Roberto Tron (tron@bu.edu)','autoTag'};
%call depPackage with file(s) from which dependences should be computed
depPackage('depPackage_test.m',opts)
