function [membershipCorrespondences]=quickshift_matchingHierarchical(X,membershipPrior,varargin)
optsMatching={};
optsScales={};
NClassesBatch=2;
%function used to compute distance
fDistance=@(x) sqrt(euclideanDistMatrix(x));

%parse optional parameters
ivarargin=1;
while ivarargin<=length(varargin)
    switch lower(varargin{ivarargin})
        case 'optsmatching'
            ivarargin=ivarargin+1;
            optsMatching=[optsMatching varargin{ivarargin}];
        case 'optsscales'
            ivarargin=ivarargin+1;
            optsScales=[optsScales varargin{ivarargin}];
        case 'nclassesbatch'
            ivarargin=ivarargin+1;
            NClassesBatch=varargin{ivarargin};
        case 'fdistance'
            ivarargin=ivarargin+1;
            fDistance=varargin{ivarargin};
        otherwise
            error(['Argument ' varargin{ivarargin} ' not valid!'])
    end
    ivarargin=ivarargin+1;
end

NPoints=size(X,2);

%divide data into batches
membershipPrior=mapValues(membershipPrior);
batchIdCnt=0;
NPriorClasses=max(membershipPrior);
batchId=NaN(size(membershipPrior));
for iClass=1:NClassesBatch:NPriorClasses
    batchIdCnt=batchIdCnt+1;
    batchIdx=iClass:min(iClass+NClassesBatch-1,NPriorClasses);
    batchId(ismember(membershipPrior,batchIdx))=batchIdCnt;
end

%stores what batches are in the first level of the hierarchy
batchIdCntFirstLevel=batchIdCnt;

%initialize structures
%points:
%   labels:             label of the cluster center they belong to
pointsLabels=1:NPoints;
%centers:
%   labels:             labels that the points use to refer to this center
%   batchId:            centers with the same batchId will be processed together
%   membershipPrior:    points with the same prior are mutually exclusive during clustering
%                       these values are computed from the batchId's of the previous level
%                       they are initialized using original membershipPrior of the points
%   scales:             at the first level of the hierarchy, they are
%   initialized to NaNs, and then the ones returned from the matching are
%   used.
centers.labels=pointsLabels;
centers.vectors=X;
centers.membershipPrior=membershipPrior;
centers.batchId=batchId;
centers.scales=NaN(size(centers.labels));

%note: batchId's and labels are incremental and not reused (even across different
%hierarchies). Moreover, batchId's at one level of are used as
%membershipPrior for the next higher level.
flagProcessNextHierarchyLevel=true;
newCentersCnt=max(centers.labels)+1;
while flagProcessNextHierarchyLevel
    batchIdList=unique(centers.batchId);
    NBatches=length(batchIdList);
    centersNew=struct('vectors',[],'labels',[],'membershipPrior',[],'scales',[]);
    for iBatch=1:NBatches
        %get batchId to be processed
        batchIdToProcess=batchIdList(iBatch);

        %see if we are at the first level of the hierarchy or not
        flagFirstLevel=batchIdToProcess<=batchIdCntFirstLevel;
        
        %get data in batch (vectors and labels)
        flagCentersBatch=centers.batchId==batchIdToProcess;
        vectorsBatch=centers.vectors(:,flagCentersBatch);
        membershipPriorBatch=centers.membershipPrior(flagCentersBatch);
        labelsBatch=centers.labels(flagCentersBatch);
        
        %if the batch contains only one class, then we do not need to run
        %the matching
        NClasses=length(unique(membershipPriorBatch));
        if NClasses==1
            %copy information to upper level
            centersNew.vectors=[centersNew.vectors vectorsBatch];
            centersNew.labels=[centersNew.labels labelsBatch];
            %compute or copy scales, depending if we are at the first level
            %of the hierarchy or not
            if flagFirstLevel
                DBatch=fDistance(vectorsBatch);
                scalesBatch=quickshift_scalesMembershipPrior(DBatch,membershipPriorBatch,optsScales{:});
            else
                scalesBatch=centers.scales(flagCentersBatch);
            end
            centersNew.scales=[centersNew.scales scalesBatch];
            flagCenters=true(1,sum(flagCentersBatch));
        else
            %compute pairwise distances
            DBatch=fDistance(vectorsBatch);

            if flagFirstLevel
                optsMatchingScales={};
            else
                optsMatchingScales={'scales',centers.scales(flagCentersBatch)};
            end
            
            %run matching on data in batch
            [labelsMatch,info]=quickshift_matching(DBatch,membershipPriorBatch,...
                optsMatching{:},'optsScales',optsScales,optsMatchingScales{:},'getcomponentdistances');

            %get vectors for new cluster centers (needs to reorder to follow
            %the labels in labelsMatch)
            flagCenters=quickshift_treeRoots(info.treeEdgesClusters);
            NCenters=sum(flagCenters);
            labelsCentersMatch=labelsMatch(flagCenters);
            newCentersVectors=vectorsBatch(:,flagCenters);
            newCentersVectors(:,labelsCentersMatch)=newCentersVectors;
            newCentersScales=info.componentDistances.default(flagCenters);
            newCentersScales(:,labelsCentersMatch)=newCentersScales;
            newCentersLabels=newCentersCnt:newCentersCnt+NCenters-1;
    
            %move forward counter for assigning labels to new centers
            newCentersCnt=newCentersCnt+NCenters;

            %collect centers in new structure 
            centersNew.vectors=[centersNew.vectors newCentersVectors];
            centersNew.labels=[centersNew.labels newCentersLabels];
            centersNew.scales=[centersNew.scales newCentersScales];

            %update data in points structure
            flagPointsBatch=ismember(pointsLabels,labelsBatch);
            pointsLabels(flagPointsBatch)=mapValues(pointsLabels(flagPointsBatch),[labelsBatch; newCentersLabels(labelsMatch)]');
        end
        centersNew.membershipPrior=[centersNew.membershipPrior batchIdToProcess*ones(1,sum(flagCenters))];
        %note: at this point the centers.batchId and points.batchId have not been updated
    end
    %break if the last number of batches was one (we reached top of hierarchy)
    if NBatches==1
        flagProcessNextHierarchyLevel=false;
    else
        
        %prepare to continue with the new level of the hierarchy        
        centers=centersNew;
        %compute new division in batches
        priorClassesList=unique(centers.membershipPrior);
        NPriorClasses=length(priorClassesList);
        batchId=NaN(size(centers.membershipPrior));
        for iClass=1:NClassesBatch:NPriorClasses
            batchIdCnt=batchIdCnt+1;
            batchIdx=priorClassesList(iClass:min(iClass+NClassesBatch-1,NPriorClasses));
            batchId(ismember(centers.membershipPrior,batchIdx))=batchIdCnt;
        end
        centers.batchId=batchId;
    end
end

membershipCorrespondences=mapValues(pointsLabels);


