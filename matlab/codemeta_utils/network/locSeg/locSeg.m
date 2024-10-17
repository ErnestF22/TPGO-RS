function locSeg(x,R,E,constraintList,constraintParametersList)
flagDisplayInfo=true;
nFigure=2;

%%
% Get dimensions of the problem
[dimData,NNodes]=size(x);
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);
NEdges=size(E,1);

%%
% Some useful definitions
uNodes=ones(NNodes,1);
INodes=eye(NNodes);

% set of indeces for getting the translation or rotation parts or both of a
% tangent vector
idxBase=(0:NNodes-1)*dimCoordinates;
idxRotation=(dimData+1:dimCoordinates)'*uNodes'+ones(dimRotations,1)*idxBase;
idxTranslation=(1:dimData)'*uNodes'+ones(dimData,1)*idxBase;
idxCoordinates=(1:dimCoordinates)'*uNodes'+ones(dimCoordinates,1)*idxBase;

%%
%

if flagDisplayInfo
    fprintf('Edges:\n')
    disp(E')
end

%%
% Obtain basis for all the shakes

J=locSegJacobianConstraintsNetwork(E,x,constraintList,constraintParametersList);
bN=null(J);
fprintfFlag(flagDisplayInfo,'Total number of shakes: %d\n', size(bN,2))


%%
% Check for common translation
HTranslation=kron(uNodes,[eye(dimData); zeros(dimRotations,dimData)]);
bNTranslationGlobal=rangeIntersection(bN, HTranslation);
fprintfFlag(flagDisplayInfo,'Rank of global translation: %d\n', size(bNTranslationGlobal,2));

%%
% Now we apply the algorithm to different centers
bNRotationGlobal=cell(1,NNodes);
bNScalingGlobal=cell(1,NNodes);
bNScalingSegments=cell(1,NNodes);
bNRotationSegments=cell(1,NNodes);

membersRotation=cell(1,NNodes);
membersScaling=cell(1,NNodes);
y=cell(1,NNodes);
for iNodeCentered=1:NNodes
    fprintfFlag(flagDisplayInfo, 'Working by centering node: %d\n', iNodeCentered);
    %compute vectors between each node and the node at the center
    y{iNodeCentered}=x-x(:,iNodeCentered)*uNodes';
    
    %%
    % Check for common rotation
    HRotation=buildCommonRotationBasis(y{iNodeCentered},R);
    bNRotationGlobal{iNodeCentered}=rangeIntersection(bN, HRotation);
    fprintfFlag(flagDisplayInfo,'\t Rank of global rotation: %d\n',...
        size(bNRotationGlobal{iNodeCentered},2));
    
    %%
    % Check for common scaling
    HScaling=reshape([y{iNodeCentered};zeros(dimRotations,NNodes)],[],1);
    bNScalingGlobal{iNodeCentered}=rangeIntersection(bN, HScaling);
    fprintfFlag(flagDisplayInfo,'\t Rank of global scaling: %d\n',...
        size(bNScalingGlobal{iNodeCentered},2));

    %%
    % Check for rotation segments
    JRotation=buildSegmentRotationConstraints(y{iNodeCentered},R,idxCoordinates);
    bNRotationSegmentsNode=null([J;JRotation]);

    tol2DRotations=1e-12;
    bNRotationSegmentsGenerated=[];
    switch dimData
        case 2
            if ~isempty(bNRotationSegmentsNode)
                C=bNRotationSegmentsNode(3:3:end,:);
                C=C';
                A=computeAllEuclideanDistancesSq(C);
                %remove row corresponding to centered node
                A(iNodeCentered,:)=[];
                newMembers=unique(A<tol2DRotations,'rows');
                NRotationSegments=size(newMembers,1);
                rankRotationSegmentsNode=size(C,1);
                CUnique=zeros(rankRotationSegmentsNode,NRotationSegments);
                for iSegment=1:NRotationSegments
                    CUnique(:,iSegment)=C(:,find(newMembers(iSegment,:),1,'first'));
                end
                %generate and check tangents
                idxMemberUnassigned=sum(newMembers,1)==0;
                %newMembers(:,iNodeCentered)=1;
%                 newMembers=locSegMarkSegmentsBoundary(newMembers,E);
                ISegments=eye(NRotationSegments);
                for iNullMode=1:rankRotationSegmentsNode
                    flagValidTangent=true(1,NRotationSegments);
                    CNullMode=CUnique(iNullMode,:);
                    for iSegment=1:NRotationSegments
                        CCentered=CNullMode-CNullMode(iSegment)*ones(1,NRotationSegments);
                        newMembers(:,idxMemberUnassigned)=ISegments(:,iSegment)*ones(1,sum(idxMemberUnassigned));
                        n=locSegBuildRotationTangentVector(y{iNodeCentered},R,CCentered,newMembers);
                        if subspace(n,bN)<1e-9
                            bNRotationSegmentsGenerated=[bNRotationSegmentsGenerated n];
                        else
                            fprintfFlag(flagDisplayInfo,'\tRejected segmentation w/ angle %s\n',num2str(subspace(n,bN)));
                            flagValidTangent(iSegment)=false;
                        end
                    end
                    if sum(flagValidTangent)<1
                        fprintfFlag(flagDisplayInfo,'\tNot enough valid tangents\n');
                        bNRotationSegmentsGenerated=[];
                        newMembers=[];
                        break
                    end
                end
                if ~isempty(newMembers)
                    newMembers(:,idxMemberUnassigned)=1;
                end
                
                newMembers(:,sum(newMembers,1)==0)=1;
                membersRotationNode{1}=newMembers;
            else
                membersRotationNode{1}=[];
            end
            bNRotationSegments{iNodeCentered}=bNRotationSegmentsGenerated;
        case 3
            fprintfFlag(flagDisplayInfo,'\t Rank of rotation segments w/ ambiguity: %d\n',...
                size(bNRotationSegmentsNode,2))

            %check for what nodes there might be an ambiguity caused by
            %global and local rotations
            bNRotationSegmentsLocal=cell(1,NNodes);
            HRotationSegmentsLocal=cell(1,NNodes);
            for iNode=1:NNodes
                HRotationSegmentsLocal{iNode}=kron(INodes(:,iNode),[zeros(3,3);eye(3)]);
                bNRotationSegmentsLocal{iNode}=rangeIntersection(bNRotationSegmentsNode,HRotationSegmentsLocal{iNode},1e-6);
            end
            
            %consider only vectors which do not contain the ambiguity
            bNRotationAmbiguity=[bNRotationSegmentsLocal{:}];
            bNRotationSegmentsNode2=rangeSubtraction(bNRotationSegmentsNode, [bNRotationAmbiguity]); %HRotation 
            rankRotationSegments2=size(bNRotationSegmentsNode2,2);
            %bNRotationSegmentsNode2=bNRotationSegmentsNode2+HRotation*randn(3,rankRotationSegments2);
            
            fprintfFlag(flagDisplayInfo,'\t Rank of rotation segments w/o ambiguity: %d\n',...
                rankRotationSegments2)

            idxENotNodeCentered=and(E(:,1)~=iNodeCentered, E(:,2)~=iNodeCentered);
            idxENodeCentered=or(E(:,1)==iNodeCentered, E(:,2)==iNodeCentered);
            ENotCentered=E(idxENotNodeCentered,:);
            ECentered=E(idxENodeCentered,:);
            membersRotationNode=cell(1,rankRotationSegments2);
            for iNullMode=1:rankRotationSegments2
                %obtain normalized infinitesimal rotation axes (in Cnorm)
                C=cnormalize(computeAllRotationIntersectionsOnEdges(bNRotationSegmentsNode2(:,iNullMode),bNRotationSegmentsLocal,y{iNodeCentered},R,ENotCentered,idxCoordinates));
                idxCValid=find(sum(C.^2)>1e-9);
                CValid=C(:,idxCValid);
                EValid=ENotCentered(idxCValid,:);
                [CUnique,membersEdgesNotCentered]=locSegUniqueNormalizedAxes(CValid);
                
                %infer clusters from grouping of edges
                NRotationSegments=size(membersEdgesNotCentered,1);
                newMembers=zeros(NRotationSegments,NNodes);
                groups=cell(1,NRotationSegments);
                for iSegment=1:NRotationSegments
                    groups{iSegment}=unique(EValid(membersEdgesNotCentered(iSegment,:)==1,:))';
                    newMembers(iSegment,groups{iSegment})=1;
                end
                
                %center can never be assigned. The current assignments are
                %from its neighbors, which can be wrong
                newMembers(:,iNodeCentered)=0;
                
%                 %infer grouping for edges connected to the center while
%                 %putting together all edge memberships
%                 membersEdge=zeros(NRotationSegments,size(E,1));
%                 membersEdge(:,idxENotNodeCentered)=membersEdgesNotCentered;
%                 for iEdgeNodeCentered=find(idxENodeCentered)'
%                     e=E(iEdgeNodeCentered,:);
%                     membersEdge(:,iEdgeNodeCentered)=newMembers(:,e(1))+newMembers(:,e(2));
%                 end
                
                %newMembers=assignCenterNode(newMembers,iNodeCentered,E);
                flagMemberUnassigned=sum(newMembers,1)==0;
                %The procedure above will not assign nodes that are only
                %connected to the center node. Get new axes from them and
                %assign them
                idxMemberUnassignedNeighbors=intersect(findNeighbors(iNodeCentered,E),find(flagMemberUnassigned));
                CNeighbors=computeRotationAxesCentered(bNRotationSegmentsNode2(:,iNullMode),y{iNodeCentered},idxMemberUnassignedNeighbors,idxTranslation);
                NUnassignedNeighbors=length(idxMemberUnassignedNeighbors);
                if NUnassignedNeighbors>0
                    newMembersUnassignedNeighbors=zeros(NUnassignedNeighbors,NNodes);
                    newMembersUnassignedNeighbors(:,idxMemberUnassignedNeighbors)=eye(NUnassignedNeighbors);
                    newMembers=[newMembers; newMembersUnassignedNeighbors];
                    CUnique=[CUnique CNeighbors];
                    groups=[groups cell(1,NUnassignedNeighbors)];
                    for iSegment=1:NUnassignedNeighbors
                        groups{iSegment+NRotationSegments}=idxMemberUnassignedNeighbors(iSegment);
                    end
                    NRotationSegments=NRotationSegments+NUnassignedNeighbors;
                    flagMemberUnassigned(idxMemberUnassignedNeighbors)=0;
                end

                %Generate tangents and check them
%                 NAxes=size(CUnique,2);
%                 for iSegment=1:NAxes
%                     groups{iSegment}=unique(EValid(membersEdgesNotCentered(iSegment,:)==1,:))';
%                 end
                [CScaled,scales,errors]=computeRotationTangentsScales(CUnique,y{iNodeCentered},bNRotationSegmentsNode2(:,iNullMode),groups);
                
                %newMembers(:,iNodeCentered)=1;
%                 newMembers=locSegMarkSegmentsBoundary(newMembers,E);
                ISegments=eye(NRotationSegments);
                flagValidTangent=true(1,NRotationSegments);
                for iSegment=1:NRotationSegments
                    CScaledCentered=CScaled-CScaled(:,iSegment)*ones(1,NRotationSegments);
                    newMembers(:,flagMemberUnassigned)=ISegments(:,iSegment)*ones(1,sum(flagMemberUnassigned));
                    n=locSegBuildRotationTangentVector(y{iNodeCentered},R,CScaledCentered,newMembers);
                    if subspace(n,bN)<1e-9
                        bNRotationSegmentsGenerated=[bNRotationSegmentsGenerated n];
                    else
                        fprintfFlag(flagDisplayInfo,'\tRejected segmentation w/ angle %s\n',num2str(subspace(n,bN)));
                        flagValidTangent(iSegment)=false;
                    end
                end
                if sum(flagValidTangent)<1
                        fprintfFlag(flagDisplayInfo,'\tNot enough valid tangents\n');
                        bNRotationSegmentsGenerated=[];
                        newMembers=[];
                end
                if ~isempty(newMembers)
                    newMembers(:,flagMemberUnassigned)=1;
                end
                membersRotationNode{iNullMode}=newMembers;
            end
            bNRotationSegments{iNodeCentered}=bNRotationSegmentsGenerated;
    end
    if ~isempty(membersRotationNode)
        if flagDisplayInfo
            fprintf('\tRotation segments by node:\n')
            cellfun(@(x) disp(x), membersRotationNode,'UniformOutput',false)
        end
    else
        fprintf('\tNo rotation segments\n')
    end
    membersRotation{iNodeCentered}=membersRotationNode;
    
    %%
    % Check for scaling segments
    JScaling=buildScalingSegmentConstraints(y{iNodeCentered},idxCoordinates);
    bNScalingSegmentsNode=null([J;JScaling]); %; bNScalingGlobal{iNodeCentered}'
    rankScalingSegmentsNode=size(bNScalingSegmentsNode,2);
    fprintfFlag(flagDisplayInfo, '\tRank of scaling segments: %d\n',rankScalingSegmentsNode)

%     if ~isempty(bNScalingSegmentsNode)
%         xdot=reshape(bNScalingSegmentsNode(idxTranslation),dimData,[]);
%         idxNodeValid=find((sum(xdot.^2))>1e-9);
%         idxEValid=and(ismember(E(:,1),idxNodeValid),ismember(E(:,2),idxNodeValid));
%         EValid=E(idxEValid,:);
%         C=POCcomputeIntersectionEdgesScaling(EValid,y{iNodeCentered},bNScalingSegmentsNode);
%         keyboard
%     end
    tolScalings=1e-8;
    bNScalingSegmentsGenerated=[];
    if ~isempty(bNScalingSegmentsNode)
        C=zeros(NNodes,size(bNScalingSegmentsNode,2));
        for iNode=1:NNodes
            if iNode~=iNodeCentered
                yi=y{iNodeCentered}(:,iNode);
                C(iNode,:)=yi'*bNScalingSegmentsNode((iNode-1)*dimCoordinates+(1:dimData),:)/(yi'*yi);
            end
        end
        C=C';
        A=computeAllEuclideanDistancesSq(C);
        %remove row corresponding to center
        A(iNodeCentered,:)=[];
        membersScalingNode=unique(A<tolScalings,'rows');
        NScalingSegments=size(membersScalingNode,1);
        CUnique=zeros(rankScalingSegmentsNode,rankScalingSegmentsNode);
        for iSegment=1:NScalingSegments
            CUnique(:,iSegment)=C(:,find(membersScalingNode(iSegment,:),1,'first'));
        end
        membersScalingNode(:,iNodeCentered)=1;
        %newMembers=assignCenterNode(newMembers,iNodeCentered,E);
        if flagDisplayInfo
            fprintf('\tScaling segments\n');
            disp(membersScalingNode)
        end
        
        %generate tangent vectors
%         for iNullMode=1:rankScalingSegmentsNode
%             n=locSegBuildScalingTangentVector(y{iNodeCentered},CUnique(iNullMode,:),membersScalingNode);
%             if subspace(n,bN)<1e-9
%                 bNScalingSegmentsGenerated=[bNScalingSegmentsGenerated n];
%             else
%                 bNScalingSegmentsGenerated=[];
%                 membersScalingNode=[];
%                 fprint('\tRejecting segmentation\n');
%                 break
%             end
%         end
%         membersScaling{iNodeCentered}=membersScalingNode;
        
        %stack all tangent vectors corresponding to translations of each
        %node. If a column is zero for all the nodes, then it means that we
        %do not have info for assigning it to a particular cluster
        xdot=reshape(permute(reshape(bNScalingSegmentsNode(idxTranslation,:),dimData,[],rankScalingSegmentsNode),[1 3 2]),rankScalingSegmentsNode*dimData,[]);%reshape(bNScalingSegmentsNode(idxTranslation),dimData,[]);
        nodeValid=find((sum(xdot.^2))>1e-9);
        
        NSegments=size(membersScalingNode,1);
        flagValidSegment=true(1,NSegments);
        for iSegment=1:NSegments
            nodeList=intersect(nodeValid,find(membersScalingNode(iSegment,:)));
            if isempty(nodeList)
                flagValidSegment(iSegment)=false;
            else
                n=POCbuildScalingTangentVector(y{iNodeCentered},nodeList);
                if subspace(n,bN)>1e-9
                    fprintfFlag(flagDisplayInfo,'\tRejected segmentation w/ angle %s\n',num2str(subspace(n,bN)));
                    flagValidSegment(iSegment)=false;
                else
                    bNScalingSegmentsGenerated=[bNScalingSegmentsGenerated n];
                end
            end
        end
        if any(~flagValidSegment)
            fprintfFlag(flagDisplayInfo,'\tRejected %d segments out of %d\n',sum(~flagValidSegment),length(flagValidSegment));
        end
        membersScaling{iNodeCentered}=membersScalingNode(flagValidSegment,:);
    end
    bNScalingSegments{iNodeCentered}=bNScalingSegmentsGenerated;
end

%Check individual rotations
bNRotationSingle=cell(NNodes,1);
I=eye(NNodes);
rankSingleRotation=zeros(1,NNodes);
for iNode=1:NNodes
    switch dimData
        case 2
            HSingleRotation=kron(I(:,iNode),[0;0;1]);
        case 3
            HSingleRotation=kron(I(:,iNode),[zeros(3); eye(3)]);
    end
    
    bNRotationSingle{iNode}=rangeIntersection(bN,HSingleRotation,1e-6);
    rankSingleRotation(iNode)=size(bNRotationSingle{iNode},2);
end
fprintf('Rank of single rotation ambiguity: %d\n',sum(rankSingleRotation))
if sum(rankSingleRotation)>0
    fprintf('    Details:\n')
    disp(rankSingleRotation)
end


%Check individual translations
switch dimData
    case 2
        HTranslationSegments=kron(eye(NNodes),[eye(2);zeros(1,2)]);
    case 3
        HTranslationSegments=kron(eye(NNodes),[eye(3);zeros(3)]);
end
bNTranslationAmbiguityFromRotations=rangeIntersection(HTranslationSegments,[bNRotationGlobal{:} bNRotationSingle{:}]);
bNTranslationSingleCheck=bN;%rangeProjection(bN,bNTranslationAmbiguityFromRotations);
bNTranslationSingle=cell(NNodes,1);
I=eye(NNodes);
rankSingleTranslation=zeros(1,NNodes);
for iNode=1:NNodes
    switch dimData
        case 2
            HSingleTraslation=kron(I(:,iNode),[eye(2); 0 0]);
        case 3
            HSingleTraslation=kron(I(:,iNode),[eye(3);zeros(3)]);
    end
    bNTranslationSingle{iNode}=rangeIntersection(bNTranslationSingleCheck,HSingleTraslation);
    rankSingleTranslation(iNode)=size(bNTranslationSingle{iNode},2);
end
fprintf('Rank of single translation ambiguity: %d\n',sum(rankSingleTranslation))
if sum(rankSingleTranslation)>0
    fprintf('    Details:\n')
    disp(rankSingleTranslation)
end


%%
% Merge different segmentations
membersRotation=flattenTwoLevelsCell(membersRotation);
if sum(cellfun(@(x) numel(x),membersRotation))>0
    [membersRotationFinal,listRotationSegments]=locSegSplitJoinSegmentations4(membersRotation);
    if flagDisplayInfo
        fprintf('Final rotation segments\n');
        cellfun(@(x) disp(x), listRotationSegments, 'UniformOutput',false)
    end
else
    membersRotationFinal=[];
end
if sum(cellfun(@(x) numel(x),membersScaling))>0
    membersScalingFinal=locSegSplitJoinSegmentations4(membersScaling);
    if flagDisplayInfo
        fprintf('Final scaling segments\n');
        disp(membersScalingFinal)
    end
else
    membersScalingFinal=[];
end
% if sum(cellfun(@(x) numel(x),membersTranslation))>0
%     membersTranslationFinal=locSegSplitJoinSegmentations4(membersTranslation);
%     if flagDisplayInfo
%         fprintf('Final translation segments\n');
%         disp(membersTranslationFinal)
%     end
% end

bNScalingSegmentsAll=[bNScalingSegments{:}];
bNRotationSegmentsAll=[bNRotationSegments{:}];

%%
% Check for translation segments
JTranslation=buildSegmentTranslationConstraints(dimData,NNodes);
bNTranslationSegments=null([J;JTranslation]);
if ~isempty(bNScalingSegmentsAll)
    switch dimData
        case 2
            HTranslationSegments=kron(membersScalingFinal',[eye(2);zeros(1,2)]);
        case 3
            HTranslationSegments=kron(membersScalingFinal',[eye(3);zeros(3)]);
    end
    bNScalingSegmentsNoTranslation=rangeSubtraction(bNScalingSegmentsAll,HTranslationSegments);
    bNScalingSegmentsTranslation=rangeIntersection(bNScalingSegmentsAll,HTranslationSegments);
else
    bNScalingSegmentsNoTranslation=[];
    bNScalingSegmentsTranslation=[];
end
if ~isempty(bNRotationSegmentsAll)
    switch dimData
        case 2
            HTranslationSegments=kron(membersRotationFinal',[eye(2);zeros(1,2)]);
        case 3
            HTranslationSegments=kron(membersRotationFinal',[eye(3);zeros(3)]);
    end
    bNRotationSegmentsNoTranslation=rangeSubtraction(bNRotationSegmentsAll,HTranslationSegments);
    bNRotationSegmentsTranslation=rangeIntersection(bNRotationSegmentsAll,HTranslationSegments);
else
    bNRotationSegmentsNoTranslation=[];
    bNRotationSegmentsTranslation=[];
end

HTranslationTruth=kron([1 1 0 0 0; 0 0 1 1 1]',[eye(3);zeros(3)]);
rankTranslationSegments=size(bNTranslationSegments,2);
fprintfFlag(flagDisplayInfo,'Rank of translation segments w/ ambiguity: %d\n',rankTranslationSegments);
bNTranslationAmbiguity=[bNScalingSegmentsAll bNScalingGlobal{:} bNRotationSegmentsAll  bNRotationSingle{:} bNRotationGlobal{:}];%      bNRotationSingle{:} bNScalingSegments{:} orth([bNRotationSegments{:}])]; %  
%bNTranslationAmbiguity=orth([bNRotationSegments{:}]);
if ~isempty(bNTranslationAmbiguity)
    bNTranslationSegments2=rangeSubtraction(bNTranslationSegments,bNTranslationAmbiguity,1e-9);
else
    bNTranslationSegments2=bNTranslationSegments;
end
bNTranslationSegments2=[bNTranslationSegments2 bNScalingSegmentsTranslation bNRotationSegmentsTranslation bNTranslationGlobal];
rankTranslationSegments2=size(bNTranslationSegments2,2);
%bNTranslationSegments2=bNTranslationSegments2+bNTranslationGlobal*randn(size(bNTranslationGlobal,2),rankTranslationSegments2);
fprintfFlag(flagDisplayInfo,'Rank of translation segments w/o ambiguity: %d\n',rankTranslationSegments2);
C=zeros(dimData*rankTranslationSegments2,NNodes);
for iNode=1:NNodes
    C(:,iNode)=reshape(bNTranslationSegments2(idxTranslation(:,iNode),:),[],1);
end
A=computeAllEuclideanDistancesSq(C);
membersTranslationFinal=unique(A<1e-9,'rows');
fprintf('Translation segments\n')
disp(membersTranslationFinal)

figure(nFigure)
nFigure=nFigure+1;
graphShow(E,x,[],'Segmentation',membersRotationFinal)
title('Rotation segmentation')

figure(nFigure)
nFigure=nFigure+1;
graphShow(E,x,[],'Segmentation',membersTranslationFinal)
title('Translation segmentation')

figure(nFigure)
nFigure=nFigure+1;
graphShow(E,x,[],'Segmentation',membersScalingFinal)
title('Scaling segmentation')

% % membersTranslation=cell(1,rankTranslationSegments2);
% for iNullMode=1:rankTranslationSegments2
%     C=cell(1,NEdges);
%     s=cell(1,NEdges);
%     for iEdge=1:NEdges
%         iNode=E(iEdge,1);
%         jNode=E(iEdge,2);
%         bNSA=[bNTranslationSegments2(:,iNullMode) bNRotationSegmentsNoTranslation];
%         bNINode=buildBaseIntersectionTranslations(iNode,bNSA,idxTranslation);
%         bNJNode=buildBaseIntersectionTranslations(jNode,bNSA,idxTranslation);
%         [C{iEdge},s{iEdge}]=rangeIntersection(bNINode,bNJNode,1e-12);
%     end
%     disp(C)
%     disp([C{:}])
% %     [~,membersEdges]=locSegUniqueNormalizedAxes(C);
% %     %infer clusters from grouping of edges
% %     NTranslationSegments=size(membersEdges,1);
% %     newMembers=zeros(NTranslationSegments,NNodes);
% %     for iSegment=1:NTranslationSegments
% %         group=unique(E(membersEdges(iSegment,:)==1,:))';
% %         newMembers(iSegment,group)=1;
% %     end
% %     membersTranslation{iNullMode}=newMembers;
% end


% %second approach
% HTranslationSegments=kron(eye(NNodes),[eye(3);zeros(3)]);
% bNTranslationAmbiguityFromRotations=rangeIntersection(HTranslationSegments,[bNRotationGlobal{:} bNRotationSingle{:}]);
% JTranslationEdges=POCBuildTranslationEdgesCheck(E);
% bNTranslationSegmentsB=null([J;JTranslation;JTranslationEdges]);
% rankTranslationSegmentsB=size(bNTranslationSegmentsB,2);
% bNTranslationSegmentsB2=rangeProjection(bNTranslationSegmentsB,[bNTranslationGlobal]); %bNTranslationAmbiguityFromRotations 
% rankTranslationSegmentsB2=size(bNTranslationSegmentsB2,2);
% 
% 

save([mfilename '_data'])
end



%%
% Auxiliary functions

%Assign the center node (which cannot be assigned with the algebraic tests)
%using the neighbors from each segment
function newMembers=assignCenterNode(members,iNodeCenter,E)
NSegments=size(members,1);
newMembers=members;
NINodeCenter=findNeighbors(iNodeCenter,E);
for iSegment=1:NSegments
    if any(members(iSegment,NINodeCenter))
        newMembers(iSegment,iNodeCenter)=1;
    else
        newMembers(iSegment,iNodeCenter)=0;
    end
end
end

% Build matrix spanning subspace of common rotations
function H=buildCommonRotationBasis(y,R)
[dimData,NNodes]=size(y);
if dimData==2
    H=[];
    S=[0 -1; 1 0];
    for iNode=1:NNodes
        H=[H; S*y(:,iNode); 1];
    end
else
    H=[];
    for iNode=1:NNodes
        H=[H; -hat(y(:,iNode)); R(:,:,iNode)'];
    end
end
end

% Build matrix for isolating tangent vectors from translations
function JCheck=buildSegmentTranslationConstraints(dimData,NNodes)
switch dimData
    case 2
        JCheck=kron(eye(NNodes),[zeros(1,2) 1]);
    case 3
        JCheck=kron(eye(NNodes),[zeros(3) eye(3)]);
end
end

% Build matrix for isolating tangent vectors from rotations
function JCheck=buildSegmentRotationConstraints(y,R,idxCoordinates)
[dimData,NNodes]=size(y);
dimConstraint=dimData;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
JCheck=zeros(dimConstraint*NNodes,dimCoordinates*NNodes);
S=[0 -1; 1 0];
idxNodes=reshape(1:NNodes*dimConstraint,dimConstraint,[]);
for iNode=1:NNodes
    idxRow=idxNodes(:,iNode);
    idxCol=idxCoordinates(:,iNode);
    switch dimData
        case 2
            JCheck(idxRow,idxCol)=[eye(2) -S*y(:,iNode)];
        case 3
            JCheck(idxRow,idxCol)=[eye(3) hat(y(:,iNode))*R(:,:,iNode)];
    end
end
end

function C=computeAllTranslationIntersectionsOnEdges(bNTranslation,E,idxTranslations)
NEdges=size(E,1);
C=cell(1,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    bNi=bNTranslation(idxTranslations(:,iNode),:);
    bNj=bNTranslation(idxTranslations(:,jNode),:);
    %vw=null([bNi bNj]);
    %C{iEdge}=vw(1:end/2,:);
    C{iEdge}=rangeIntersection(bNi,bNj);
end
end

function C=computeRotationAxesCentered(nSeg,y,neighbors,idxTranslation)
NNeighbors=length(neighbors);
dimData=size(y,1);
if dimData~=3
    error('This function is not implemented for 2D data')
end
C=zeros(dimData,NNeighbors);
for iiNode=1:NNeighbors
    iNode=neighbors(iiNode);
    C(:,iiNode)=cnormalize(cross(y(:,iNode),nSeg(idxTranslation(:,iNode))));
end
end
    
function C=computeAllRotationIntersectionsOnEdges(nSeg,bNAmbiguity,y,R,E,idxCoordinates)
NEdges=size(E,1);
dimData=3;

C=zeros(3,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    idx=[idxCoordinates(:,iNode); idxCoordinates(:,jNode)];
    
    H=[-hat(y(:,iNode)); R(:,:,iNode)';-hat(y(:,jNode)); R(:,:,jNode)'];
    bNij=[nSeg(idx) bNAmbiguity{iNode}(idx,:) bNAmbiguity{jNode}(idx,:)];

    vw=null([H -bNij]);
    
    c=vw(1:3,:);
    if ~isempty(c)
        C(:,iEdge)=c;
    else
        C(:,iEdge)=zeros(dimData,1);
    end
end
end

function JCheck=buildScalingSegmentConstraints(y,idxCoordinates)
[dimData,NNodes]=size(y);
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);
dimConstraint=dimCoordinates-1;
normV=sqrt(sum(y.^2));
flagZeroNormV=normV<1e-12;

%for those nodes which have y{iNodeCentered} with zero norm, we will need constraints of
%dimension dimConstraint+1 (i.e., dimCoordinates)
JCheck=zeros(dimConstraint*NNodes+sum(flagZeroNormV),dimCoordinates*NNodes);

idxRowJCheck=0;
for iNode=1:NNodes
    idx=idxCoordinates(:,iNode);
    if flagZeroNormV(iNode)
        %effectively freeze center node
        JCheck(idxRowJCheck+(1:dimCoordinates),idx)=eye(dimCoordinates);
        idxRowJCheck=idxRowJCheck+dimCoordinates;
    else
        %build orthogonal complement of y(:,iNode) using Householder
        %tranformations
        orthConstr=computeOrthogonalComplement(y(:,iNode));

        JCheck(idxRowJCheck+(1:dimConstraint),idx)=...
            blkdiag(orthConstr',eye(dimRotations));
        idxRowJCheck=idxRowJCheck+dimConstraint;
    end
end
end

function bNINode=buildBaseIntersectionTranslations(iNode,bNSA,idxTranslation)
idxINode=idxTranslation(:,iNode);
bNINode=[bNSA(idxINode,:)];%[bNTranslationSegments2(idxINode,iNullMode) bNScalingSegments(idxINode,:) bNGlobalScaling(idxINode,:)];
end

function [vScaled,scales,errors]=computeRotationTangentsScales(vNorm,y,n,groups)
NNodes=size(y,2);
NSegments=length(groups);
scales=zeros(1,NSegments);
vScaled=zeros(size(vNorm));
dimData=3;
dimCoordinates=6;
idxTranslation=(1:dimData)'*ones(1,NNodes)+ones(dimData,1)*(0:NNodes-1)*dimCoordinates;

for iSegment=1:NSegments
    A=[];
    b=[];
    for g=groups{iSegment}
        A=[A; -hat(y(:,g))*vNorm(:,iSegment)];
        b=[b; n(idxTranslation(:,g))];
    end
    scales(iSegment)=(A'*b)/(A'*A);
    vScaled(:,iSegment)=scales(iSegment)*vNorm(:,iSegment);
    if nargout>2
        errors(iSegment)=norm(A*scales(iSegment)-b);
    end
end
end