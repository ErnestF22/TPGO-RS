%Test idea for a factorization-based full pose bearing-only localization
function POCFactorizationBearingLocalization
flagDirectFix=false;
flagSoftFix=true;
flagScaleInequalities=false;
flagCheckIdxMatrix=true;
flagCheckMeasurements=true;
flagCheckBearingConstraints=true;
flagCheckRotationColumn=true;

methodPerturbation='small';

nbNodes=5;

A=adjgallery(nbNodes,'kneigh',2);
t_node=testNetworkBuildTestNetwork('A',A);

E=testNetworkGetEdges(t_node);
Ri=t_node.Ri;
Ti=t_node.Ti;
tij=cnormalize(t_node.Tij);

%change coordinates to fix one of the nodes with pose (I,0)
iEdgeFix=4;
iNodeFix=E(iEdgeFix,1);
jNodeFix=E(iEdgeFix,2);

[Ri,Ti]=RTFix(Ri,Ti,iNodeFix,'references');

WInfo=lowRankLocalization_infoInit(nbNodes,'rotationAugmented','inodefix',iNodeFix,'dim',3);

W=lowRankLocalization_groundTruthW(Ri,Ti,WInfo);

idxMatrix=WInfo.idxMatrix;

if flagCheckIdxMatrix
    disp('Check of makeIdxMatrix')
    iEdge=2;
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    TjNode=Ti(:,jNode);
    RiNode=Ri(:,:,iNode);
    disp(' - Check entry for one edge endpoint')
    disp([RiNode'*TjNode W(idxMatrix(:,iNode,jNode))])
    disp(' - Check that top right block is the identity')
    disp(W(squeeze(idxMatrix(:,WInfo.iNodeFix,WInfo.jNodeR))))
    disp(' - Check that fourth eval is zero')
    s=svd(W);
    disp(s(1:4)')
end

if flagCheckBearingConstraints
    disp('Check fundamental relation between directions of the measurements and W')
    iEdge=2;
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    tijEdge=tij(:,iEdge);
    TiNode=Ti(:,iNode);
    TjNode=Ti(:,jNode);
    RiNode=Ri(:,:,iNode);
    %disp([cnormalize(Tj-Ti) Ri*tij])
    %disp(orthComplementProjector(Ri*tij)*(Tj-Ti))
    %disp(orthComplementProjector(tij)*Ri'*(Tj-Ti))
    disp(' - Bearing measurements')
    disp(orthComplement(tijEdge)'*...
        (W(idxMatrix(:,iNode,jNode))-W(idxMatrix(:,iNode,iNode))))
end
if flagCheckRotationColumn
    disp(' - Rotations in the last block column')
    e=lowRankLocalization_errorOrthonormalityRotations(WInfo,W);
    disp(e)
end

if false
    disp('Check of fundamental relation between scale of the measurements and W')
    iEdge=4;
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    tijEdge=tij(:,iEdge);
    TiNode=Ti(:,iNode);
    TjNode=Ti(:,jNode);
    RiNode=Ri(:,:,iNode);
    %disp([cnormalize(Tj-Ti) Ri*tij])
    %disp(orthComplementProjector(Ri*tij)*(Tj-Ti))
    %disp(orthComplementProjector(tij)*Ri'*(Tj-Ti))
    disp(tijEdge'*(RiNode'*TjNode-RiNode'*TiNode))
    
    aij=makeSingleConstraintScale(tijEdge,iNode,jNode,nbNodes);
    disp(aij*W(:))
end

%matrix from measurements
[AMeasurements,bMeasurements]=lowRankLocalization_bearings_constraints(WInfo,tij,E);

if flagSoftFix
    AFixTranslation=zeros(3*nbNodes,3*nbNodes*(nbNodes+3));
    AFixTranslation(sub2ind(size(AFixTranslation),(1:3*nbNodes)',vec(idxMatrix(:,:,iNodeFix))))=1;
    bFixTranslation=zeros(3*nbNodes,1);
    AFixRotation=zeros(9,3*nbNodes*(nbNodes+3));
    AFixRotation(sub2ind(size(AFixRotation),(1:9)',vec(idxMatrix(:,iNodeFix,WInfo.jNodeR))))=1;
    bFixRotation=vec(eye(3));
    AFixScale=zeros(1,3*nbNodes*(nbNodes+3));
    AFixScale(:,idxMatrix(:,iNodeFix,jNodeFix))=tij(:,iEdgeFix)';
    AFixScale(:,idxMatrix(:,iNodeFix,iNodeFix))=-tij(:,iEdgeFix)';
    bFixScale=norm(Ti(:,jNodeFix)-Ti(:,iNodeFix));
    
    AMeasurements=cat(1,AMeasurements,AFixTranslation,AFixRotation,AFixScale);
    bMeasurements=cat(1,bMeasurements,bFixTranslation,bFixRotation,bFixScale);
end

if flagScaleInequalities
    AScales=zeros(nbEdges,3*nbNodes*(nbNodes+3));
    for iEdge=1:nbEdges
        iNode=E(iEdge,1);
        jNode=E(iEdge,2);
        tijEdge=tij(:,iEdge);
        aij=makeSingleConstraintScale(tijEdge,iNode,jNode,nbNodes);
        AScales(iEdge,:)=aij;
        bScales=11.3137*ones(nbEdges,1);
    end
else
    AScales=[];
    bScales=[];
end

if flagDirectFix
    AFix=zeros(6,3*nbNodes*(nbNodes+3));
    bFix=zeros(6,1);
    AFix(1:3,idxMatrix(:,iNodeFix,iNodeFix))=eye(3);
    bFix(1:3)=W(idxMatrix(:,iNodeFix,iNodeFix));
    AFix(4:6,idxMatrix(:,iNodeFix,jNodeFix))=eye(3);
    bFix(4:6)=W(idxMatrix(:,iNodeFix,jNodeFix));

    %Remove edge used to fix ambiguity from measurements
    AMeasurements(idxA(:,iEdgeFix),:)=[];
    bMeasurements(idxA(:,iEdgeFix))=[];
    
    AMeasurements=[AMeasurements; AFix];
    bMeasurements=[bMeasurements; bFix];
end

if flagCheckMeasurements
    disp('Error from measurement constraints')
    disp(norm(AMeasurements*W(:)-bMeasurements,1))
end

%add perturbation to initial condition
switch methodPerturbation
    case 'none'
        WInitial=W;
    case 'subspace'
        WInitial=W+10*reshape(AMeasurements'*randn(nbConstraints,1),sz(1),sz(2));
    case 'full'
        WInitial=W+10*randn(size(W));
    case 'small'
        WInitial=W+0.1*randn(size(W));
    case 'empty'
        %Ask for automatic initialization
        WInitial=[];
    otherwise
        error('methodPerturbation not recognized')
end        

solver='admm';

switch lower(solver)
    case 'cvx'
        cvx_begin
            variable WEstimated(3*nbNodes,nbNodes+3)
            f=norm_nuc(WEstimated);
            minimize(f)
            subject to
                AMeasurements*vec(WEstimated)==bMeasurements
                %consequence of fixing one of the nodes at the origin
                WEstimated(vec(idxMatrix(:,:,iNodeFix)))==0
                %consequence of fixing one of the rotations to the identity
                WEstimated(squeeze(idxMatrix(:,iNodeFix,jNodeR)))==eye(3)
                %scale constraints
                AScales*WEstimated(:)>=bScales
        %         %very simple relaxation of the convex hull of SO(3)
        %         for iNode=1:nbNodes
        %             Wi=WEstimated(squeeze(idxMatrix(:,iNode,jNodeR)));
        %             Wi<=1
        %             Wi>=-1
        %             sum(abs(vec(Wi)))<=5
        %         end
        cvx_end
    case 'admm'
        [WEstimated,~,output]=fixedRankLSQP(AMeasurements,bMeasurements,...
            AScales,bScales,size(W),3,...
            'maxIt',2000,'initialSolution',WInitial,...
            'referenceSolution',W);
        figure(1)
        fixedRankLSQP_plotErrors(output)
end
disp([W WEstimated])
[RRef,TRef]=recoverFactorization(W);
[REstimated,TEstimated]=recoverFactorization(WEstimated);
save([mfilename '_data'])
keyboard

function [R,T]=recoverFactorization(W)
[U,S,V]=svd(W,'econ');
R=matUnstack(U(:,1:3));
T=S(1:3,1:3)*V(:,1:3)';
alignment=T(:,end-2:end);
R=matUnstack(matStack(R)*alignment);
T=alignment\T;

function aij=makeSingleConstraintScale(tij,iNode,jNode,nbNodes)
idxMatrix=makeIdxMatrix(nbNodes);
aij=zeros(1,3*nbNodes*(nbNodes+3));
aij(:,idxMatrix(:,iNode,iNode))=-tij';
aij(:,idxMatrix(:,iNode,jNode))=tij';
