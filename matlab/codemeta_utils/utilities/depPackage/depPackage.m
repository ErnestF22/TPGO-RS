%function depPackage(filename,opts)
%Computes dependencies, stage them, add autorights, add additional files,
%and create zip file
function depPackage(filename,opts)
flagStop=false;

%defaults
if ~isfield(opts,'autorightsPreambleFile')
    opts.autorightsPreambleFile=fullfile(mfilepath(),'licensePreamble.txt');
end
if ~isfield(opts,'additionalFiles')
    opts.additionalFiles={};
end  
if ~any(contains(lower(opts.additionalFiles),'license'))
    %add a license file if not present (according to the criteria above)
    opts.additionalFiles=[opts.additionalFiles fullfile(mfilepath(),'LICENSE')];
end
    

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
        %normalize dir name
        opts.baseDir=realpath(opts.baseDir);
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
