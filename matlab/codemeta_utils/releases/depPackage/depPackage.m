%function depPackage(filename,opts)
%Computes dependencies, stage them, add autorights, add additional files,
%and create zip file
function depPackage(filename,opts)

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

% Ver: 09-Jul-2018 12:54:16

% Roberto Tron (tron@bu.edu)

flagStop=false;

%defaults
if ~isfield(opts,'autorightsPreambleFile')
    opts.autorightsPreambleFile=fullfile(mfilepath(),'licensePreamble.txt');
end
if ~isfield(opts,'additionalFiles') || ~any(contains(lower(opts.additionalFiles),'license'))
    %add a license file if not present (according to the criteria above)
    opts.additionalFiles={fullfile(mfilepath(),'LICENSE')};
end
    

%normalize dir name
opts.baseDir=realpath(opts.baseDir);

%check if directory exists and ask user what to do
if exist(opts.baseDir,'dir')
    a='';
    while length(a)~=1 || (a~='a' && a~='n' && a~='r')
        a=lower(input(['Directory ' opts.baseDir ' exists. [A]dd and proceed/[R]emove and proceed/Do [N]ot proceed. '],'s'));
    end
    if a=='n'
        flagStop=true;
    end
    if a=='r'
        cmd=['find "' opts.baseDir '" -not -path ''*/\.*'' -delete'];
        system(cmd);
    end
end

if ~flagStop
    [oldNameList,newNameList]=depPackageStaging(filename,opts.pairs,opts.baseDir);
    depPackageAddAutorights(newNameList,opts.autorightsPreambleFile,opts.autorightsOpts{:});

    if isfield(opts,'additionalFiles')
        disp('Copying additional files')
        copyAdditionalFiles(opts)
    end
    
    disp('Creating ZIP archive')
    cdOld=cd(opts.baseDir);
    cmd=['zip -r ' opts.packageName ' *'];
    [b,output]=system(cmd);
    if b>1
        disp(output)
        error('Error creating zip archive')
    end

    cd(cdOld)

    disp(['Archive ' fullfile(opts.baseDir,opts.packageName) '.zip created'])
end

function copyAdditionalFiles(opts)
for iFile=1:length(opts.additionalFiles)
    copyfile(opts.additionalFiles{iFile},opts.baseDir);
end
