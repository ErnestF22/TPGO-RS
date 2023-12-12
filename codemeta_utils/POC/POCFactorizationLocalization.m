function POCFactorizationLocalization
resetRands()
flagCheckTranslationMeasurements=false;
sigmaNoise=0.0;

nbNodes=5;
A=adjgallery(nbNodes,'kneigh',2);
t_node=testNetworkBuildTestNetwork('A',A);

E=testNetworkGetEdges(t_node);

Ri=t_node.Ri;
Ti=t_node.Ti;
Tij=t_node.Tij.*repmat(t_node.lambdaij,3,1);
Tij=Tij+sigmaNoise*randn(size(Tij));

%change coordinates to fix one of the nodes with pose (I,0)
iNodeFix=1;
[Ri,Ti]=RTFix(Ri,Ti,iNodeFix,'references');

t_node.Ritruth=Ri;
t_node.Titruth=Ti;
WInfo=lowRankLocalization_infoInit(nbNodes,'rotationAugmented','inodefix',iNodeFix,'dim',3);
W=lowRankLocalization_groundTruthW(Ri,Ti,WInfo);

idxMatrix=lowRankLocalization_idxMatrix(WInfo);


if flagCheckTranslationMeasurements
    iEdge=4;
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    RiNode=Ri(:,:,iNode);
    TiNode=Ti(:,iNode);
    TjNode=Ti(:,jNode);
    disp('Check fundamental relation of translation measurements')
    disp('- Check definition relative translation measurements')
    disp(Tij(:,iEdge)-RiNode'*(TjNode-TiNode))
    disp('- Check with W matrix')
    disp(Tij(:,iEdge)-(W(idxMatrix(:,iNode,jNode))-W(idxMatrix(:,iNode,iNode))))
end

[ATranslations,bTranslations]=lowRankLocalization_translations_constraints(WInfo,Tij,E);
disp(norm(ATranslations*vec(W)-bTranslations))

[AFix,bFix]=lowRankLocalization_gauge_constraints(WInfo);
disp(norm(AFix*vec(W)-bFix))

AMeasurements=[ATranslations;AFix];
bMeasurements=[bTranslations;bFix];

WInitial=[];

[WEstimated,~,output]=fixedRankLSQP(AMeasurements,bMeasurements,...
    [],[],size(W),3,...
    'maxIt',300,'initialSolution',WInitial,...
    'referenceSolution',W);
%figure(1)
%fixedRankLSQP_plotErrors(output)

[RiEstimated,TiEstimated]=lowRankLocalization_solution_extractProjection(WEstimated,WInfo);
disp('Error RiEstimated-Ri')
disp(norm(vec(RiEstimated-Ri)))
disp('Error TiEstimated-Ti')
disp(norm(vec(TiEstimated-Ti)))

save([mfilename '_data'])
%keyboard