function lowRankLocalization_gauge_constraints_test
nbNodes=5;
A=adjgallery(nbNodes,'kneigh',2);
t_node=testNetworkBuildTestNetwork('A',A);

E=testNetworkGetEdges(t_node);
Ri=t_node.Ri;
Ti=t_node.Ti;


%change coordinates to fix one of the nodes with pose (I,0)
iNodeFix=1;
[Ri,Ti]=RTFix(Ri,Ti,iNodeFix,'references');

WInfo=lowRankLocalization_infoInit(nbNodes,'rotationAugmented','inodefix',iNodeFix,'dim',3);
W=lowRankLocalization_groundTruthW(Ri,Ti,WInfo);

[AFix,bFix]=lowRankLocalization_gauge_constraints(WInfo);
disp(norm(AFix*vec(W)-bFix))
