%function depPackageStaging(filename,pairs,baseDir)
%Takes the matlab file FILENAME, runs depPackageNameList and then copies
%all the files to BASEDIR
function [oldNameList,newNameList]=depPackageStaging(filename,pairs,baseDir)

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

% Ver: 04-Jul-2018 13:12:59

% Roberto Tron (tron@bu.edu)


disp('Generating dependencies')
[newNameList,oldNameList]=depPackageNameList(filename,pairs);
disp('Copying files')
for iName=1:length(newNameList)
    destDir=fullfile(baseDir,fileparts(newNameList{iName}));
    if ~exist(destDir,'dir')
        b=mkdir(destDir);
        if ~b
            error('Error creating directory %s for staging')
        end
    end
    %fprintf('\t%s\n',newNameList{iName})
    newNameList{iName}=fullfile(baseDir,newNameList{iName});
    b=copyfile(oldNameList{iName},newNameList{iName});
    if ~b
        error('Error staging the file %s',newNameList{iName})
    end
end
