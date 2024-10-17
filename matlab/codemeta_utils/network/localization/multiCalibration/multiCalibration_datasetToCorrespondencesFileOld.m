%function c=multiCalbration_datasetToCorrespondencesFile(dirName)
%Looks for files in the given directory containing a pair of the following
%variables: apriltag_planes, velodyne_planes, lgrass_planes, rgrass_planes,
%hokuyo_endpts. For each pair found, it creates a struct element in the
%array c which contains the file name, variable names, contained data and
%type of the correspondence (plane-plane or plane-point).
%
%If present, this function also loads the file nodePairs.txt, which
%contains triplets of <fileName> <nodeName1> <nodeName2>, which indicate,
%for each file in the dataset, the human-readable name of the nodes
%involved in the correspondences. Note that for plane-point correspondence,
%the name of the node producing points is always assumed to come second.
%This information is added to the struct element of the cell array.
%
function c=multiCalibration_datasetToCorrespondencesFile(dirName)
flagNodePairs=false;
fileNameNodePairs=fullfile(dirName,'nodePairs.txt');
if exist(fileNameNodePairs)
    fid=fopen(fileNameNodePairs,'r');
    cNodePairs=textscan(fid,'%s %s %s');
    fclose(fid);
    flagNodePairs=true;
end

d=dir(fullfile(dirName,'*.mat'));
NFiles=length(d);
c={};
for iFiles=1:NFiles
    detectedFieldType={};
    detectedFieldName={};
    fileName=fullfile(dirName,d(iFiles).name);
    s=load(fileName);
    if isfield(s,'apriltag_planes')
        detectedFieldType{end+1}='plane';
        detectedFieldName{end+1}='apriltag_planes';
    end
    if isfield(s,'velodyne_planes')
        detectedFieldType{end+1}='plane';
        detectedFieldName{end+1}='velodyne_planes';
    end
    if isfield(s,'lgrass_planes')
        detectedFieldType{end+1}='plane';
        detectedFieldName{end+1}='lgrass_planes';
    end
    if isfield(s,'rgrass_planes')
        detectedFieldType{end+1}='plane';
        detectedFieldName{end+1}='rgrass_planes';
    end
    if isfield(s,'hokuyo_endpts')
        detectedFieldType{end+1}='point';
        detectedFieldName{end+1}='hokuyo_endpts';
    end

    if length(detectedFieldName)==2
        %Assume fields with points are detected after those with plane
        cCurrent=[];
        cCurrent.fileName=fileName;
        if flagNodePairs
            [flagMember,idxMember]=ismember(d(iFiles).name,cNodePairs{1});
            if flagMember
                cCurrent.nodeNames={cNodePairs{2}{idxMember} cNodePairs{3}{idxMember}};
            end
        end
        if ~strcmp(detectedFieldType{2},'point')
            cCurrent.type='plane-plane';
        else
            cCurrent.type='plane-point';
        end
        cCurrent.fieldNames=detectedFieldName;
        cCurrent.rawData{1}=s.(detectedFieldName{1});
        cCurrent.rawData{2}=s.(detectedFieldName{2});
        
        c{end+1}=cCurrent;
    else
        disp(['Skipped ' fileName])
    end
end

fprintf('Number of files: %d\n',NFiles)
fprintf('Number of correspondences: %d\n',length(c))

