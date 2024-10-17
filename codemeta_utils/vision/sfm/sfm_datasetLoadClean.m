%Load a dataset from a data file and with noiseless projected points
%function data=sfm_datasetLoadClean(datasetName)
function data=sfm_datasetLoadClean(datasetName)
switch datasetName
    case 'synthetic'
        data=sfm_datasetGenerate();
    otherwise
        thresholdNCorrespondences=15;
        methodRebuildMatches='structureID';
        
        s=load(['sfmdata_' datasetName]);
        dataOriginal=s.data;
        
        %copy over some fields
        data.imgIdName=dataOriginal.imgIdName;
        data.imgSize=dataOriginal.imgSize;
        data.resizeFactor=dataOriginal.resizeFactor;
        data.imgFileName=dataOriginal.imgFileName;
        data.poseTruth=dataOriginal.poseTruth;
        data.projection=dataOriginal.projection;
        data.calibration=dataOriginal.calibration;
        
        switch lower(methodRebuildMatches)
            case 'reprojections'
                [feature,match]=rebuildMatchesUsingReprojections(dataOriginal);
            case 'structureid'
                [feature,match]=rebuildMatchesUsingStructureIDs(dataOriginal);
            otherwise
                error('methodRebuildMatches not valid')
        end
        
        data.feature=feature;
        data=sfm_featureNormalize(data);

        %remove matches that now do not have enough correspondences
        flagMatch=cellfun(@(x) size(x,2),{match.idxMatch})>thresholdNCorrespondences;
        match=match(flagMatch);
        
        data.match=match;
end
data=sfm_matchPoseTruth(data,'memberMatch','match');

function [feature,match]=rebuildMatchesUsingReprojections(data)
%associate each feature to the closest reprojection of a 3-D point which is
%less than a threshold of pixels away
threshold=2;
NFeatures=length(data.feature);
feature=repmat(struct('location',[]),1,NFeatures);
idxValid=cell(NFeatures,1);
idxStructure=cell(NFeatures,1);

X=data.structureFiltered.location;
for iFeature=1:NFeatures
    P=data.projection(:,:,iFeature);
    xReproj=projectFromP(P,X);
    x=data.feature(iFeature).location;
    d=euclideanDistMatrix(xReproj,x);
    [dMin,idxMin]=min(d);
    flagValid=dMin<threshold;
    feature(iFeature).location=xReproj(:,idxMin(flagValid));
    
    %keep structure membership and indeces of valid points for later
    idxStructure{iFeature}=idxMin;
    idxValid{iFeature}=find(flagValid);
end

%modify the matches to reflect the changes in the available features
NMatches=length(data.match);
match=repmat(struct('idxImg',[],'idxMatch',[]),...
    1,NMatches);
for iMatch=1:NMatches
    thisMatch=data.match(iMatch);
    idxImg=thisMatch.idxImg;
    idxMatch=thisMatch.idxMatch;
    
    %remove matches that point to features without or with different 3-D points
    flagValidMatch=true(1,size(idxMatch,2));
    idxMatchStructure=zeros(size(idxMatch));
    for k=1:2
        idxk=idxValid{idxImg(k)};
        flagValidMatch=flagValidMatch & ismember(idxMatch(k,:),idxk);
        idxMatchStructure(k,:)=idxStructure{idxImg(k)}(idxMatch(k,:));
    end
    flagValidMatch=flagValidMatch & (idxMatchStructure(1,:)==idxMatchStructure(2,:));
    
    idxMatch=idxMatch(:,flagValidMatch);
        
    match(iMatch).idxImg=idxImg;
    match(iMatch).idxMatch=zeros(size(idxMatch));

    %loop across the two images
    for k=1:2
        idxk=idxValid{idxImg(k)};
        match(iMatch).idxMatch(k,:)=mapValues(idxMatch(k,:),idxk');
    end
end

function [feature,match]=rebuildMatchesUsingStructureIDs(data)
%filter features based on whether they have a 3-D point or not, and
%store the corresponding index information for later
NFeatures=length(data.feature);
feature=repmat(struct('location',[],'locationNormalized',[]),...
    1,NFeatures);
idxValid=cell(NFeatures,1);
for iFeature=1:NFeatures
    idx=find(data.feature(iFeature).structureFilteredMembership>0);
    feature(iFeature).location=data.feature(iFeature).locationReprojected(:,idx);
    feature(iFeature).locationNormalized=data.feature(iFeature).locationReprojectedNormalized(:,idx);
    idxValid{iFeature}=idx;
end

%modify the matches to reflect the changes in the available
%features
NMatches=length(data.match);
match=repmat(struct('idxImg',[],'idxMatch',[]),...
    1,NMatches);
for iMatch=1:NMatches
    thisMatch=data.match(iMatch);
    m=thisMatch.idxMatch;
    idxImg=thisMatch.idxImg;

    %get stucture membership for features referenced by the matches
    structId=zeros(size(m));
    for k=1:2
        structId(k,:)=data.feature(idxImg(k)).structureFilteredMembership(m(k,:));
    end

    %keep only matches that correspond to features of valid 3-D
    %points with the same structure membership
    flagThisMatch=all(structId>0) & structId(1,:)==structId(2,:);
    thisMatch=sfm_rawFilterDataWithFlag(thisMatch,flagThisMatch);

    match(iMatch).idxImg=thisMatch.idxImg;
    match(iMatch).idxMatch=zeros(size(thisMatch.idxMatch));

    for k=1:2
        idxk=idxValid{idxImg(k)};
        match(iMatch).idxMatch(k,:)=mapValues(thisMatch.idxMatch(k,:),idxk');
    end
end

