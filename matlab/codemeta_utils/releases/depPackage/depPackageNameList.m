%function [oldNameList,newNameList]=depPackageNameList(filename,pairs,baseDir)
%Takes the matlab file FILENAME, runs FDEP and runs DEPPACKAGESUBSTITUTIONS
%using pairs.
%If FILENAME is a cell array of strings, the dependencies are obtained by
%running FDEP on each entry and then merging the results.
function [newNameList,oldNameList]=depPackageNameList(filename,pairs)

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

% Ver: 04-Jul-2018 13:12:03

% Roberto Tron (tron@bu.edu)

if ~iscell(filename)
    filename={filename};
end

oldNameList={};
for iName=1:length(filename)
    %%This was required pre-MATLAB2015a
    %[filePath,fileNameNoExt]=fileparts(filename{iName});
    %cdold=cd;
    %if ~isempty(filePath)
    %    if ~exist(filePath,'dir')
    %        error('File name path "%s" not found',filePath); 
    %    end
    %    cd(filePath);
    %end
    %rp=fdep(fileNameNoExt,'-q');   
    %cd(cdold);
    %oldNameList=[oldNameList;rp.fun];
    oldNameList=[oldNameList matlab.codetools.requiredFilesAndProducts(filename{iName})'];
end
oldNameList=unique(oldNameList);
newNameList=depPackageSubstitutions(oldNameList,pairs,'removeRoot','ensureFileSep');
