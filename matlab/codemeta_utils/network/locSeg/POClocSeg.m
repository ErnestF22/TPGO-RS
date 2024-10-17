function POClocSeg
resetRands()
dimData=3;
graphSampleNum=23;
constraintSampleNum=3;
nodeCentered=3;

[E,x,R,constraintList,constraintParametersList]=locSegSampleNetworkGallery(graphSampleNum,'dimData',dimData,'constraintNum',constraintSampleNum);
figure(1)
graphShow(E,x,R)

NNodes=max(E(:));
NEdges=size(E,1);
if NNodes<3
    nodeCentered=1;
end
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);


% set of indeces for getting the translation or rotation parts or both of a
% tangent vector
uNodes=ones(NNodes,1);
INodes=eye(NNodes);

idxBase=(0:NNodes-1)*dimCoordinates;
idxRotation=(dimData+1:dimCoordinates)'*uNodes'+ones(dimRotations,1)*idxBase;
idxTranslation=(1:dimData)'*uNodes'+ones(dimData,1)*idxBase;
idxCoordinates=(1:dimCoordinates)'*uNodes'+ones(dimCoordinates,1)*idxBase;


J=locSegJacobianConstraintsNetwork(E,x,constraintList,constraintParametersList);

%cell array of cell arrays containing inputs for locSegApplyTransformation
transformationList={};


%x=x-x(:,nodeCentered)*ones(1,NNodes);

switch dimData
    case 2
        dimCoordinates=3;
    case 3
        dimCoordinates=6;
end

%Check for global translations
bN=null(J);
fprintf('Total number of shakes: %d\n', size(bN,2))
fprintf('Rank of J: %d\n', rank(J))
fprintf('Size of J:')
disp(size(J))

HCommonTranslation=buildMatrixCommonTranslationCheck(dimData,dimCoordinates,NNodes);
% wv=null([bN -HCommonTranslation]);
% bNGlobalTranslation=bN*wv(1:size(bN,2),:);
bNGlobalTranslation=rangeIntersection(bN, HCommonTranslation,1e-6);
[dxGlobalTranslation,dRGlobalTranslation]=locSegNullspace2dxdR(bNGlobalTranslation,R,dimData,dimCoordinates);
rankGlobalTranslation=size(bNGlobalTranslation,2);

for iRank=1:rankGlobalTranslation
    TTransformation=bNGlobalTranslation(1:dimData,iRank);
    transformationList{end+1}={'translation',{TTransformation},1:NNodes};
end

fprintf('Rank of global translation ambiguity: %d\n', rankGlobalTranslation);

%From now on, we consider transformations centered at one of the nodes
y=centerVector(x,nodeCentered);
%JCenter=buildConstraintsCenter(dimData,NNodes,nodeCentered);

%Check for global rotation
HRotation=buildHatMatrixCommonRotationCheck(y,R);
% wv=null([bN HRotation]);
% bNGlobalRotation=bN*wv(1:size(bN,2),:);
bNGlobalRotation=rangeIntersection(bN, HRotation,1e-6);
rankGlobalRotation=size(bNGlobalRotation,2);
[dxCommonRotation,dRCommonRotation]=locSegNullspace2dxdR(bNGlobalRotation,R,dimData,dimCoordinates);

cTransformation=x(:,nodeCentered);
for iRank=1:rankGlobalRotation
    vTransformation=R(:,:,1)*bNGlobalRotation(dimData+1:dimCoordinates,iRank);
    transformationList{end+1}={'rotation',{vTransformation, cTransformation},1:NNodes};
end

fprintf('Rank of global rotation ambiguity: %d\n', rankGlobalRotation);

%Check for global scalings
%JFixRotation=buildConstraintsFixedRotation(dimData,NNodes);
%bNGlobalScaling=null([J;JCenter;JFixRotation]);
%wv=null([bNGlobalScaling -HGlobalScaling]);
%bNGlobalScaling=bNGlobalScaling*wv(1:end-1,:);
HGlobalScaling=reshape([y;zeros(dimCoordinates-dimData,NNodes)],[],1);
% wv=null([bN -HGlobalScaling]);
% bNGlobalScaling=bN*wv(1:end-1,:);
bNGlobalScaling=rangeIntersection(bN, HGlobalScaling);
rankGlobalScaling=size(bNGlobalScaling,2);
[dxGlobalScaling,dRGlobalScaling]=locSegNullspace2dxdR(bNGlobalScaling,R,dimData,dimCoordinates);

if rankGlobalScaling>0
    cTransformation=x(:,nodeCentered);
    aTransformation=1;
    transformationList{end+1}={'scaling',{aTransformation, cTransformation},1:NNodes};
end

fprintf('Rank of global scaling ambiguity: %d\n', rankGlobalScaling)

%Check for scaling segments
JScalingCheck=buildConstraintsScaleCheck(y);
bNScalingSegments=null([J;JScalingCheck;bNGlobalScaling']);
[dxScalingSegments,dRScalingSegments]=locSegNullspace2dxdR(bNScalingSegments,R,dimData,dimCoordinates);

fprintf('Rank of cluster scaling ambiguity: %d\n',size(bNScalingSegments,2))

tolScalings=1e-8;
if ~isempty(bNScalingSegments)
    C=zeros(NNodes,size(bNScalingSegments,2));
    for iNode=1:NNodes
        if iNode~=nodeCentered
            yi=y(:,iNode);
            C(iNode,:)=yi'*bNScalingSegments((iNode-1)*dimCoordinates+(1:dimData),:)/(yi'*yi);
        end
    end
    A=computeAllEuclideanDistancesSq(C');
    %remove row corresponding to center
    A(nodeCentered,:)=[];
    membersScaling=unique(A<tolScalings,'rows');
    NScalingSegments=size(membersScaling,1);
    fprintf('Number of scaling clusters: %d\n',NScalingSegments)
    disp(membersScaling)
    
    cTransformation=x(:,nodeCentered);
    aTransformation=1;
    for iSegment=1:NScalingSegments
        nodeListSegment=find(membersScaling(iSegment,1));
        transformationList{end+1}={'scale',{aTransformation, cTransformation}, nodeListSegment};
    end
    
%     for iNullMode=1:size(bNScalingSegments,2)
%         fprintf('    Clustering according to non-global scalings, mode %d\n',iNullMode);
%         n=dxScalingSegments(:,:,iNullMode);
%         gammas=sum(n.*n)./sum(n.*y);
%         A=abs(gammas([1:nodeCentered-1 nodeCentered+1:NNodes])'*ones(1,NNodes)...
%             -ones(NNodes-1,1)*gammas);
%         disp(A<tolGamma)
%     end
        
        
%     idxTranslations=(1:dimData)'*ones(1,NNodes)+dimCoordinates*ones(dimData,1)*(0:NNodes-1);
%     tolAngle=1e-6;
% 
%     for iDim=1:dimData
%         A=computeAllPairSemiAngles(bNScalingSegments(idxTranslations(iDim,:),:)');
%         fprintf('Clustering according to scaling, dim. %d\n',iDim)
%         disp(A([1:nodeCentered-1 nodeCentered+1:NNodes],:)<tolAngle)
%     end
end

%Check rotation segments
JRotationCheck=buildConstraintsRotationCheck2(y,R);
%JRotationCheckB=buildConstraintsRotationCheck2b(y,R);
%JGlobalRotation=buildConstraintsRotationCheck(y);
%bNRotationSegments=null([J;JCenter;JGlobalRotation;JFixRotation;bNGlobalRotation']);
bNRotationSegments=null([J;JRotationCheck]);
[dxRotationSegments,dRRotationSegments]=locSegNullspace2dxdR(bNRotationSegments,R,dimData,dimCoordinates);
idxENotNodeCentered=and(E(:,1)~=nodeCentered, E(:,2)~=nodeCentered);

fprintf('Rank of cluster rotation ambiguity: %d\n',size(bNRotationSegments,2))

tolLocalRotationsSubspaceAngle=1e-12;
tol2DRotations=1e-12;
tol3DRotations=1e-12;
switch dimData
    case 2
%             gammas=sum(y.*([0 -1;1 0]*dxRotationSegments(:,:,iNullMode)))./(sum(y.*y));
%             A=abs(gammas([1:nodeCentered-1 nodeCentered+1:NNodes])'*ones(1,NNodes)...
%                 -ones(NNodes-1,1)*gammas);
        C=bNRotationSegments(3:3:end,:);
        A=computeAllEuclideanDistancesSq(C');
        disp(A<tol2DRotations)
    case 3
        %check for what nodes there might be an ambiguity caused by local
        %rotations
        ENotCentered=E(idxENotNodeCentered,:);
        
        for iNullMode=1:size(bNRotationSegments,2)
            C3{iNullMode}=cnormalize(computeAllCrossProductsOnEdges(dxRotationSegments(:,:,iNullMode),ENotCentered));
            [C3norm{iNullMode},flagsC3norm]=locSegUniqueNormalizedAxes(C3{iNullMode});
        end
        
        I=eye(NNodes);
        vLocalRotations=zeros(size(y));
        bNRotationSegmentsLocal=cell(1,NNodes);
        HRotationSegmentsLocal=cell(1,NNodes);
        for iNode=1:NNodes
            HRotationSegmentsLocal{iNode}=kron(I(:,iNode),[zeros(3,3);eye(3)]);
            bNRotationSegmentsLocal{iNode}=rangeIntersection(bNRotationSegments,HRotationSegmentsLocal{iNode});
        end
        rankRotationSegmentsLocal=cellfun(@(x) size(x,2), bNRotationSegmentsLocal);
            
        a=zeros(1,NNodes);
        for iNode=1:NNodes
            if iNode~=nodeCentered
                a(iNode)=subspace(bNRotationSegments,kron(I(:,iNode),[zeros(3,1);R(:,:,iNode)'*y(:,iNode)]));
                if a(iNode)<tolLocalRotationsSubspaceAngle
                    vLocalRotations(:,iNode)=y(:,iNode);
                end
            end
        end
        bNRotationAmbiguityWG=[HRotation bNRotationSegmentsLocal{:}];
        %bNSegWG=rangeSubtraction(bNRotationSegments, bNRotationAmbiguityWG,1e-6);
        bNRotationAmbiguity=[bNRotationSegmentsLocal{:}];
        bNSeg=rangeSubtraction(bNRotationSegments, [HRotation bNRotationAmbiguity],1e-6);
        groups={};
        for iNullMode=1:size(bNSeg,2)
            [dxSeg,dRSeg]=locSegNullspace2dxdR(bNSeg(:,iNullMode),R,dimData,dimCoordinates);
            C2prime{iNullMode}=cnormalize(computeAllCrossProductsOnEdges(dxSeg,ENotCentered));
            C2{iNullMode}=cnormalize(computeAllIntersectionsOnEdges2(bNSeg(:,iNullMode),bNRotationSegmentsLocal,y,R,ENotCentered));
            [C2norm{iNullMode},flagsC2norm]=locSegUniqueNormalizedAxes(C2{iNullMode});
            disp([E(idxENotNodeCentered,:)'; C2{iNullMode}]);
            for iSegment=1:size(C2norm{iNullMode},2)
                groups{iSegment}=unique(ENotCentered(flagsC2norm(iSegment,:)==1,:))';
            end
            [C2Scaled{iNullMode},scales,errors]=computeRotationTangentsScales(C2norm{iNullMode},y,bNSeg(:,iNullMode),groups);
            A2{iNullMode}=computeTranslationBackError(C2Scaled{iNullMode},y,bNSeg(:,iNullMode));
        end
        JEqualRot=[];
        for iSegment=1:length(groups)
            g=[nodeCentered groups{iSegment}];
            %nEqualRot=generateTangentVectorFromGroupRotationConstraints(bNSeg,bNRotationSegmentsLocal,g,R);
            nEqualRot{iSegment}=generateTangentVectorFromGroupRotationConstraintsFix(bNSeg,bNRotationSegmentsLocal,HRotation,g,R);
        end
        fprintf('Number of rotation clusters (latest method): %d\n',length(groups))
%         EE=extractTransitiveEdgePairs(E);
%         C3=computeAllIntersectionsOnTransitiveEdges(bNSeg(:,1),bNRotationAmbiguity,y,R,EE);
%         displayValidAxes(EE,C3)
        
%         for iNullMode=1:size(bNRotationSegments,2)
%             %C=computeAllCrossProductsOnEdges(dxRotationSegments(:,:,iNullMode),E(idxENotNodeCentered,:));
%             w=extractAxesFromRotationTangent(bNRotationSegments(:,iNullMode),R);
%             C=computeAllIntersectionsOnEdges(vLocalRotations,w,ENotCentered);
%             Cnorm=locSegUniqueNormalizedAxes(C);
%             A1=abs(Cnorm'*dxRotationSegments(:,:,iNullMode));
%             A2=checkRotationIntersectionCompatibility(Cnorm,vLocalRotations,w);
%             newMembers=and(A1<tol3DRotations,A2<tol3DRotations);
%             newMembers=assignCenterNode(newMembers,nodeCentered,E);
%             newAttributes=cellfun(@(x) {x}, mat2cell(Cnorm,size(Cnorm,1),ones(size(Cnorm,2),1)),'UniformOutput',false);
%             if iNullMode==1
%                 membersRotationSegments=newMembers;
%                 attributesRotationSegments=newAttributes;
%             else
%                 [membersRotationSegments,attributesRotationSegments]=locSegSplitJoinSegmentations(membersRotationSegments,newMembers,attributesRotationSegments,newAttributes);
%             end
%             
% %             fprintf('    Clustering according to rotations, null mode %d\n',iNullMode)
% %             disp(newMembers)
%             
%         end
end
% fprintf('Number of rotation clusters: %d\n',size(membersRotationSegments,1))
% disp(membersRotationSegments)
% disp(attributesRotationSegments)

%Check individual rotations
bNSingleRotation=cell(NNodes,1);
I=eye(NNodes);
rankSingleRotation=zeros(1,NNodes);
for iNode=1:NNodes
    switch dimData
        case 2
            HSingleRotation=kron(I(:,iNode),[0;0;1]);
        case 3
            HSingleRotation=kron(I(:,iNode),[zeros(3); eye(3)]);
    end
    
    bNSingleRotation{iNode}=rangeIntersection(bN,HSingleRotation,1e-6);
    rankSingleRotation(iNode)=size(bNSingleRotation{iNode},2);
end
fprintf('Rank of single rotation ambiguity: %d\n',sum(rankSingleRotation))
if sum(rankSingleRotation)>0
    fprintf('    Details:\n')
    disp(rankSingleRotation)
end


%Check individual translations
bNSingleTranslation=cell(NNodes,1);
I=eye(NNodes);
rankSingleTranslation=zeros(1,NNodes);
JTranslation=buildSegmentTranslationConstraints(dimData,NNodes);
bNSingleTranslationCheck=null([J;JTranslation]);
for iNode=1:NNodes
    switch dimData
        case 2
            HSingleTraslation=kron(I(:,iNode),[eye(2); 0 0]);
        case 3
            HSingleTraslation=kron(I(:,iNode),[eye(3);zeros(3)]);
    end
    bNSingleTranslation{iNode}=rangeIntersection(bN,HSingleTraslation);
    %bNSingleTranslation{iNode}=rangeIntersection(rangeSubtraction(bNSingleTranslationCheck,bNSingleRotation{iNode}),HSingleTraslation);
    rankSingleTranslation(iNode)=size(bNSingleTranslation{iNode},2);
end
fprintf('Rank of single translation ambiguity: %d\n',sum(rankSingleTranslation))
if sum(rankSingleTranslation)>0
    fprintf('    Details:\n')
    disp(rankSingleTranslation)
end

% JTranslation=buildSegmentTranslationConstraints(dimData,NNodes);
% bNTranslationSegments=null([J;JTranslation]);%;bNGlobalTranslation';[bNSingleTranslation{:}]';bNScalingSegments']);
% bNTranslationAmbiguity=[bNGlobalTranslation bNGlobalScaling bNScalingSegments bNSingleTranslation{:}]; %
% bNTranslationSegments2=rangeSubtraction(bNTranslationSegments,bNTranslationAmbiguity);
% rankTranslationSegments=size(bNTranslationSegments,2);
% rankTranslationSegments2=size(bNTranslationSegments2,2);
% bNTranslationSegments2=bNTranslationSegments2+bNGlobalTranslation*randn(size(bNGlobalTranslation,2),rankTranslationSegments2);
% fprintf('Rank of translation segments w/ ambiguity: %d\n',rankTranslationSegments);
% fprintf('Rank of translation segments w/o ambiguity: %d\n',rankTranslationSegments2);
% 
% if rankTranslationSegments2>0
%     C=cell(1,NEdges);
%     s=cell(1,NEdges);
%     for iEdge=1:NEdges
%         iNullMode=1;
%         iNode=E(iEdge,1);
%         jNode=E(iEdge,2);
%         bNSA=[bNTranslationSegments2(:,iNullMode) bNScalingSegments];
%         bNINode=buildBaseIntersectionTranslations(iNode,bNSA,bNSingleTranslation,idxTranslation);
%         bNJNode=buildBaseIntersectionTranslations(jNode,bNSA,bNSingleTranslation,idxTranslation);
% 
%         [C{iEdge},s{iEdge}]=rangeIntersection(bNINode,bNJNode,1e-12);
%     end
%     disp(E')
%     disp(C)
%     disp([C{:}])
% end
JTranslation=buildSegmentTranslationConstraints(dimData,NNodes);
bNTranslationSegments=null([J;JTranslation]);%;bNGlobalTranslation';[bNSingleTranslation{:}]';bNScalingSegments']);
%bNTranslationSegments=[bNSingleTranslation{:}];
bNTranslationAmbiguity=[bNGlobalTranslation bNGlobalRotation bNGlobalScaling bNScalingSegments bNSeg bNSingleRotation{:}];%
bNTranslationSegments2=rangeSubtraction(bNTranslationSegments,bNTranslationAmbiguity);
rankTranslationSegments=size(bNTranslationSegments,2);
rankTranslationSegments2=size(bNTranslationSegments2,2);
bNTranslationSegments2=bNTranslationSegments2+bNGlobalTranslation*randn(size(bNGlobalTranslation,2),rankTranslationSegments2);
fprintf('Rank of translation segments w/ ambiguity: %d\n',rankTranslationSegments);
fprintf('Rank of translation segments w/o ambiguity: %d\n',rankTranslationSegments2);

if rankTranslationSegments2>0
    C=cell(1,NEdges);
    s=cell(1,NEdges);
    for iEdge=1:NEdges
        iNullMode=2;
        iNode=E(iEdge,1);
        jNode=E(iEdge,2);
        bNSA=[bNTranslationSegments2(:,iNullMode) bNScalingSegments];
        bNINode=buildBaseIntersectionTranslations(iNode,bNSA,bNSingleTranslation,idxTranslation);
        bNJNode=buildBaseIntersectionTranslations(jNode,bNSA,bNSingleTranslation,idxTranslation);

        [C{iEdge},s{iEdge}]=rangeIntersection(bNINode,bNJNode,1e-12);
    end
    disp(E')
    disp(C)
    disp([C{:}])
end
%disp([E';C{:}])
%keyboard
%return
% M=zeros(dimData*rankTranslationSegments2*NNodes,rankTranslationSegments2+size(bNTranslationAmbiguity,2));
% idxM=reshape(1:size(M,1),dimData*rankTranslationSegments2,NNodes);
% for iNode=1:NNodes
%     idxINode=idxTranslation(:,iNode);
%     MTemp=[];
%     for iNullMode=1:rankTranslationSegments2
%         MTemp=blkdiag(MTemp,bNTranslationSegments2(idxINode,iNullMode));
%     end
%     M(idxM(:,iNode),:)=[MTemp,...
%         kron(ones(rankTranslationSegments2,1),bNTranslationAmbiguity(idxINode,:))];
% end
% 
% C=cell(1,NEdges);
% for iEdge=1:NEdges
%     iNode=E(iEdge,1);
%     jNode=E(iEdge,2);
%     C{iEdge}=rangeIntersection(M(idxM(:,iNode),:),M(idxM(:,jNode)));
% end

tolTranslations=1e-8;
C=zeros(dimData*rankTranslationSegments,NNodes);
for iNode=1:NNodes
    C(:,iNode)=reshape(bNTranslationSegments(idxTranslation(:,iNode),:),[],1);
end
A=computeAllEuclideanDistancesSq(C);
membersTranslation=unique(A<tolTranslations,'rows');
fprintf('Segmentation according to translations:\n')
disp(membersTranslation)



%bNDetectedNoTranslationSegments=[bNGlobalTranslation bNGlobalRotation bNGlobalScaling bNScalingSegments bNSeg bNSingleRotation{:}];
bNDetected=[bNGlobalTranslation bNTranslationSegments2 bNGlobalRotation bNGlobalScaling bNScalingSegments bNSeg bNSingleRotation{:}];% bNSingleTranslation{:}];

bNResidual=rangeSubtraction(bN,bNDetected,1e-6);

fprintf('Subspace distance between bN and bNDetected: %.4e\n',subspace(bN,bNDetected))
fprintf('Number of redundant modes: %d\n', size(null(bNDetected),2))
fprintf('Number of unknown modes: %d\n',size(bNResidual,2))

%[U,S,V]=svd(bNDetected,'econ');
% bNResidual=bN-bNDetected*(bNDetected\bN);
% disp('Rank of residual shakes')
% disp(sum(svd(bNResidual)>1e-8))
%keyboard

% nodeCentered2=4;
% nodeCentered2=2;

% bNCentered=locSegCenterBasis(bN,dimCoordinates,dimData,nodeCentered);
% bNCenteredNormalized=cnormalize(bNCentered')';
% [dx,dR]=nullspace2dxdR(bNCentered,R,dimData,dimCoordinates);
% disp('bNCenteredNormalized')
% displayMatrixInRowBlocks(bNCenteredNormalized,dimCoordinates)

% bNCentered2=locSegBiCenterBasis(bN,dimCoordinates,dimData,nodeCentered,nodeCentered2);
% bNCentered2Normalized=cnormalize(bNCentered2')';
% [dx2,dR2]=nullspace2dxdR(bNCentered2,R,dimData,dimCoordinates);
% 
% disp('bNCenteredNormalized2')
% displayMatrixInRowBlocks(bNCentered2Normalized,dimCoordinates)

% if dimData==3 && (graphSampleNum==8 || graphSampleNum==10) && constraintSampleNum==2
%     bV=reshape([y;zeros(3,NNodes)],[],1);
%     bN=null([J;bV']);
%     bNCentered2=locSegBiCenterBasis(bN,dimCoordinates,dimData,3,1);
%     [dxCentered2,dRCentered2]=nullspace2dxdR(bNCentered2,R,dimData,dimCoordinates);
%     %[dx,dR]=nullspace2dxdR(bN,dimData,R,dimCoordinates);
%     if size(bNCentered2,2)>0
%         nullMode=1;graphAnimate(E,x,dxCentered2(:,:,nullMode),R,dRCentered2(:,:,nullMode))
%     else
%         disp('No twisting modes')
%     end
%     %keyboard
% end

save([mfilename '_data'])

function A=computeTranslationBackError(vScaled,y,n)
NNodes=size(y,2);
NSegments=size(vScaled,2);
dimData=3;
dimCoordinates=6;
idxTranslation=(1:dimData)'*ones(1,NNodes)+ones(dimData,1)*(0:NNodes-1)*dimCoordinates;
A=zeros(NSegments,NNodes);

for iNode=1:NNodes
    for iSegment=1:NSegments
        A(iSegment,iNode)=norm(n(idxTranslation(:,iNode))+hat(y(:,iNode))*vScaled(:,iSegment));
    end
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


function members=assignCenterNode(members,iCenterNode,E)
[NSegments,NNodes]=size(members);
NICenterNode=findNeighbors(iCenterNode,E);
for iSegment=1:NSegments
    if any(members(iSegment,NICenterNode))
        members(iSegment,iCenterNode)=1;
    else
        members(iSegment,iCenterNode)=0;
    end
end

function H=buildMatrixCommonTranslationCheck(dimData,dimCoordinates,NNodes)
H1=[eye(dimData); zeros(dimCoordinates-dimData,dimData)];
H=kron(ones(NNodes,1), H1);

function A=computeAllPairSemiAngles(y)
tolNorm=1e-12;
NVectors=size(y,2);
A=zeros(NVectors);
for iVector=1:NVectors
    for jVector=1:NVectors
        yi=y(:,iVector);
        yj=y(:,jVector);
        nvi=norm(yi);
        nvj=norm(yj);
        if nvi<tolNorm || nvj<tolNorm
            A(iVector,jVector)=0;
        else
            A(iVector,jVector)=acos(min(1,abs((yi/nvi)'*(yj/nvj))));
        end
    end
end

function C=computeAllCrossProductsOnEdges(y,E)
%Computes all the cross products between all the pairs specified
%in E of vectors in y
NEdges=size(E,1);
C=zeros(3,NEdges);
for iEdge=1:NEdges
    yi=y(:,E(iEdge,1));
    yj=y(:,E(iEdge,2));
    
    C(:,iEdge)=cross(yi,yj);
end

function C=computeAllIntersectionsOnEdges(y,w,E)
%Computes all the cross products between all the pairs specified
%in E of vectors in y
NEdges=size(E,1);
dimData=size(y,1);
C=zeros(3,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    yi=y(:,iNode);
    yj=y(:,jNode);
    wi=w(:,iNode);
    wj=w(:,jNode);
    
    c=rangeIntersection([yi wi],[yj wj],1e-6);
    if ~isempty(c)
        C(:,iEdge)=c;
    else
        C(:,iEdge)=zeros(dimData,1);
    end
end

function C=computeAllIntersectionsOnEdges2(nSeg,bNAmbiguity,y,R,E)
NEdges=size(E,1);
NNodes=size(R,3);
dimData=3;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
%indeces for each node in a tangent vector
idxBN=ones(dimCoordinates,1)*((0:NNodes-1)*dimCoordinates)+(1:dimCoordinates)'*ones(1,NNodes);

% %preprocess bNAmbiguity and nSeg to remove the effect of each rotation from
% %the rotation part of the tangent vectors
% for iNode=1:NNodes
%     idx=idxBN(:,iNode);
%     bNAmbiguity(idx,:)=R(:,:,iNode)*bNAmbiguity(idx,:);
%     nSeg(idx,:)=R(:,:,iNode)*nSeg(idx,:);
% end

C=zeros(3,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    idx=[idxBN(:,iNode); idxBN(:,jNode)];
    
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

function nEqualRot=generateTangentVectorFromGroupRotationConstraints(bNSeg,bNRotationSegmentsLocal,g,R)
dimData=size(R,1);
NNodes=size(R,3);
dimCoordinates=locSegDimData2DimCoordinates(dimData);

idxRot=(dimData+1:dimCoordinates)'*ones(1,NNodes)+ones(dimCoordinates-dimData,1)*((0:NNodes-1)*dimCoordinates);
JEqualRot=[];
for idx=1:length(g)-1
    JEqualRotRow=zeros(3,NNodes*dimCoordinates);
    iNode=g(idx);
    jNode=g(idx+1);
    JEqualRotRow(:,idxRot(:,iNode))=R(:,:,iNode);
    JEqualRotRow(:,idxRot(:,jNode))=-R(:,:,jNode);
    JEqualRot=[JEqualRot; JEqualRotRow];
end
nnEqualRot=null(JEqualRot*[bNSeg bNRotationSegmentsLocal{g}]);
nEqualRot=[bNSeg bNRotationSegmentsLocal{g}]*nnEqualRot;

function nEqualRot=generateTangentVectorFromGroupRotationConstraintsFix(bNSeg,bNRotationSegmentsLocal,HRotation,g,R)
dimData=size(R,1);
NNodes=size(R,3);
dimCoordinates=locSegDimData2DimCoordinates(dimData);

idxRot=(dimData+1:dimCoordinates)'*ones(1,NNodes)+ones(dimCoordinates-dimData,1)*((0:NNodes-1)*dimCoordinates);
JEqualRot=[];
for iNode=g
    JEqualRotRow=zeros(3,NNodes*dimCoordinates);
    JEqualRotRow(:,idxRot(:,iNode))=R(:,:,iNode);
    JEqualRot=[JEqualRot; JEqualRotRow];
end
bNAllowedRotation=[bNSeg bNRotationSegmentsLocal{g} HRotation];
nnEqualRot=null(JEqualRot*bNAllowedRotation);
nEqualRot=bNAllowedRotation*nnEqualRot;

function C=computeAllIntersectionsOnTransitiveEdges(nSeg,bNAmbiguity,y,R,EE)
NEdges=size(EE,1);
NNodes=size(R,3);
dimData=3;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
%indeces for each node in a tangent vector
idxBN=ones(dimCoordinates,1)*((0:NNodes-1)*dimCoordinates)+(1:dimCoordinates)'*ones(1,NNodes);

C=cell(1,NEdges);
for iEdge=1:NEdges
    iNode=EE(iEdge,1);
    jNode=EE(iEdge,2);
    kNode=EE(iEdge,3);
    idx=[idxBN(:,iNode); idxBN(:,jNode); idxBN(:,kNode)];
    
    H=[-hat(y(:,iNode)); R(:,:,iNode)';-hat(y(:,jNode)); R(:,:,jNode)';-hat(y(:,kNode)); R(:,:,kNode)'];
    bNij=orth([nSeg(idx) bNAmbiguity(idx,:)]);

    vw=null([H -bNij]);
    
    C{iEdge}=cnormalize(vw(1:3,:));
end

function displayValidAxes(EE,C)
NEdges=size(EE,1);
D=[];
for iEdge=1:NEdges
    if ~isempty(C{iEdge})
        D=[D [EE(iEdge,:)'; C{iEdge}]];
    end
end
disp(D)
    
function y=centerVector(x,nodeCentered)
y=x-x(:,nodeCentered)*ones(1,size(x,2));

function JCheck=buildConstraintsScaleCheck(y)
[dimData,NNodes]=size(y);
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);
dimConstraint=dimCoordinates-1;
normV=sqrt(sum(y.^2));
flagZeroNormV=normV<1e-12;
%for those nodes which have y with zero norm, we will need constraints of
%dimension dimConstraint+1 (i.e., dimCoordinates)
JCheck=zeros(dimConstraint*NNodes+sum(flagZeroNormV),dimCoordinates*NNodes);
idxRowJCheck=0;
for iNode=1:NNodes
    if flagZeroNormV(iNode)
        %effectively freeze center node
        JCheck(idxRowJCheck+(1:dimCoordinates),dimCoordinates*(iNode-1)+(1:dimCoordinates))=eye(dimCoordinates);
        idxRowJCheck=idxRowJCheck+dimCoordinates;
    else
        %build orthogonal complement of y(:,iNode) using Householder
        %tranformations
        
        %compute d=y/norm(y)-eN
        d=cnormalize(y(:,iNode));
        d(end)=d(end)-1;
        d=cnormalize(d);
        orthConstr=eye(dimData,dimData-1)-2*d*d(1:end-1)';
        %orthConstr=null(y(:,iNode)');

        JCheck(idxRowJCheck+(1:dimConstraint),dimCoordinates*(iNode-1)+(1:dimCoordinates))=blkdiag(orthConstr',eye(dimRotations));
        idxRowJCheck=idxRowJCheck+dimConstraint;
    end
end

function JCheck=buildConstraintsRotationCheck(y)
[dimData,NNodes]=size(y);
dimConstraint=1;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
JCheck=zeros(dimConstraint*NNodes,dimCoordinates*NNodes);
for iNode=1:NNodes
    if norm(y(:,iNode))>1e-12
        orthConstr=y(:,iNode)';
        JCheck(dimConstraint*(iNode-1)+(1:dimConstraint),dimCoordinates*(iNode-1)+(1:dimData))=orthConstr';
    end
end


function w=extractAxesFromRotationTangent(n,R)
dimData=size(R,1);
[dimCoordinates,dimRotations]=locSegDimData2DimCoordinates(dimData);
NNodes=size(n,1)/dimCoordinates;
w=zeros(dimRotations,NNodes);
switch dimData
    case 2
        error('Not implemented yet!')
    case 3
        for iNode=1:NNodes
            w(:,iNode)=R(:,:,iNode)*n(dimCoordinates*(iNode-1)+(dimData+1:dimCoordinates),:);
        end
end

function n=generateGlobalRotationTangent(y,R,w)
if ~exist('w','var')
    w=cnormalize(randn(3,1));
end
[dimData,NNodes]=size(y);
dimCoordinates=locSegDimData2DimCoordinates(dimData);
n=[];
for iNode=1:NNodes
    n=[n; -hat(y(:,iNode))*w; R(:,:,iNode)'*w];
end
    
function JCheck=buildConstraintsRotationCheck2(y,R)
[dimData,NNodes]=size(y);
dimConstraint=dimData;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
JCheck=zeros(dimConstraint*NNodes,dimCoordinates*NNodes);
S=[0 -1; 1 0];
for iNode=1:NNodes
    switch dimData
        case 2
            JCheck(dimConstraint*(iNode-1)+(1:dimConstraint),dimCoordinates*(iNode-1)+(1:dimCoordinates))=[eye(2) -S*y(:,iNode)];
        case 3
            JCheck(dimConstraint*(iNode-1)+(1:dimConstraint),dimCoordinates*(iNode-1)+(1:dimCoordinates))=[eye(3) hat(y(:,iNode))*R(:,:,iNode)];
    end
end

function JCheck=buildConstraintsRotationCheck2b(y,R)
[dimData,NNodes]=size(y);
dimConstraint=1;
dimCoordinates=locSegDimData2DimCoordinates(dimData);
JCheck=zeros(dimConstraint*NNodes,dimCoordinates*NNodes);
S=[0 -1; 1 0];
for iNode=1:NNodes
    switch dimData
        case 2
            error('Not implemented yet')
            %JCheck(dimConstraint*(iNode-1)+(1:dimConstraint),dimCoordinates*(iNode-1)+(1:dimCoordinates))=[eye(2) -S*y(:,iNode)];
        case 3
            JCheck(dimConstraint*(iNode-1)+(1:dimConstraint),dimCoordinates*(iNode-1)+(dimData+1:dimCoordinates))=y(:,iNode)'*R(:,:,iNode);
    end
end

function H=buildHatMatrixCommonRotationCheck(y,R)
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


function JCenter=buildConstraintsCenter(dimData,NNodes,nodeCentered)
dimCoordinates=locSegDimData2DimCoordinates(dimData);
JCenter=zeros(dimData,dimCoordinates*NNodes);
JCenter(:,dimCoordinates*(nodeCentered-1)+(1:dimData))=eye(dimData);

function JFixRotation=buildConstraintsFixedRotation(dimData,NNodes)
[~,dimRotations]=locSegDimData2DimCoordinates(dimData);
JFixRotation=kron(eye(NNodes),[zeros(dimRotations,dimData) eye(dimRotations)]);

function A=checkRotationIntersectionCompatibility(Cnorm,y,w)
NK=size(Cnorm,2);
NNodes=size(y,2);
A=zeros(NK,NNodes);
for iK=1:NK
    for iN=1:NNodes
        A(iK,iN)=subspace(Cnorm(:,iK),[y(:,iN) w(:,iN)]);
    end
end

function JCheck=buildSegmentTranslationConstraints(dimData,NNodes)
switch dimData
    case 2
        JCheck=kron(eye(NNodes),[zeros(2,1) 1]);
    case 3
        JCheck=kron(eye(NNodes),[zeros(3) eye(3)]);
end

function C=computeAllTranslationIntersectionsOnEdges(bNTranslation,E,idxTranslations)
NEdges=size(E,1);
C=cell(1,NEdges);
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    bNi=bNTranslation(idxTranslations(:,iNode),:);
    bNj=bNTranslation(idxTranslations(:,jNode),:);
    vw=null([bNi bNj]);
    C{iEdge}=vw(1:end/2,:);
end

function bNINode=buildBaseIntersectionTranslations(iNode,bNSA,bNSingleTranslation,idxTranslation)
idxINode=idxTranslation(:,iNode);
bNINode=[bNSA(idxINode,:) bNSingleTranslation{iNode}(idxINode,:)];%[bNTranslationSegments2(idxINode,iNullMode) bNScalingSegments(idxINode,:) bNGlobalScaling(idxINode,:)];
% if ~isempty(bNSingleTranslation{iNode})
%     bNINode=[bNINode bNSingleTranslation{iNode}(idxINode,:)];
% end        
