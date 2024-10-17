%Traverse feature matching graph and give feature membership for 3-D points
%function data=sfm_structureExtractMembership(data)
%Requires fields feature().matchMembership*. Adds field structure with
%sub-fields feature and idxFeature. Uses
%feature().matchMembershipFlagVisited for doing the graph sweep
function data=sfm_structureExtractMembership(data)
NFeatures=length(data.feature);

%allocate space for visit flags
for iFeature=1:NFeatures
    NPoints=size(data.feature(iFeature).location,2);
    NList=data.feature(iFeature).matchMembershipMaxCount;
    data.feature(iFeature).matchMembershipFlagVisited=false(NList,NPoints);
    %put to true all the out-of-list indicators
    for iPoint=1:NPoints
        data.feature(iFeature).matchMembershipFlagVisited(data.feature(iFeature).matchMembershipCount(iPoint)+1:end,iPoint)=true;
    end
    data.feature(iFeature).structureMembership=zeros(1,NPoints);
end

X=[];
idxX=1;
%do graph sweep
%loop over images
for iFeature=1:NFeatures
    %find possible roots, i.e., features that appear in matches that have
    %not been accounted for
    NPoints=size(data.feature(iFeature).location,2);
    for iRoot=1:NPoints
        if ~all(data.feature(iFeature).matchMembershipFlagVisited(:,iRoot))
            %init lists of feature groups and indeces relative to this root
            %(will grow during the graph sweep and then stored in the
            %structure)
            rootFeature=iFeature;
            rootIdxFeature=iRoot;

            %init visit queue
            rootNextFlag=~data.feature(iFeature).matchMembershipFlagVisited(:,iRoot);
            rootNextFeature=data.feature(iFeature).matchMembershipListFeature(rootNextFlag,iRoot);
            rootNextIdxFeature=data.feature(iFeature).matchMembershipListIdxFeature(rootNextFlag,iRoot);

            %visit children
            while ~isempty(rootNextFeature)

                %get front of visit queue
                childFeature=rootNextFeature(1);
                childIdxFeature=rootNextIdxFeature(1);

                %check which entries of the child correspond to features
                %already visited or that are scheduled to be visited
                childFlagVisited=data.feature(childFeature).matchMembershipFlagVisited(:,childIdxFeature) ...
                    | (ismember(data.feature(childFeature).matchMembershipListFeature(:,childIdxFeature), [rootFeature;rootNextFeature]) ...
                    & ismember(data.feature(childFeature).matchMembershipListIdxFeature(:,childIdxFeature), [rootIdxFeature;rootNextIdxFeature]));

                %mark entries in child
                data.feature(childFeature).matchMembershipFlagVisited(childFlagVisited,childIdxFeature)=true;

                %add whatever remains to the visit queue
                rootNextFeature=[rootNextFeature; ...
                    data.feature(childFeature).matchMembershipListFeature(~childFlagVisited,childIdxFeature)]; %#ok<AGROW>
                rootNextIdxFeature=[rootNextIdxFeature; ...
                    data.feature(childFeature).matchMembershipListIdxFeature(~childFlagVisited,childIdxFeature)]; %#ok<AGROW>

                %add child to root's lists
                rootFeature=[rootFeature; childFeature]; %#ok<AGROW>
                rootIdxFeature=[rootIdxFeature; childIdxFeature]; %#ok<AGROW>

                %remove head of visit queue
                rootNextFeature(1)=[];
                rootNextIdxFeature(1)=[];
            end
            
            %TODO: right now the flags at the root are not marked
            
            
            %add feature to the structure array
            X(idxX).feature=rootFeature; %#ok<AGROW>
            X(idxX).idxFeature=rootIdxFeature;
            idxX=length(X)+1;
        end
    end
end

%compact the structure into matrices
featureCount=arrayfun(@(x) length(x.feature),X);
featureMaxCount=max(featureCount);
NPoints=length(X);

structure=struct('feature',zeros(featureMaxCount,NPoints),...
    'idxFeature',zeros(featureMaxCount,NPoints),...
    'featureCount',featureCount,...
    'featureMaxCount',featureMaxCount);

for iPoint=1:NPoints
    nFeature=featureCount(iPoint);
    structure.feature(1:nFeature,iPoint)=X(iPoint).feature;
    structure.idxFeature(1:nFeature,iPoint)=X(iPoint).idxFeature;
end

data.structure=structure;
