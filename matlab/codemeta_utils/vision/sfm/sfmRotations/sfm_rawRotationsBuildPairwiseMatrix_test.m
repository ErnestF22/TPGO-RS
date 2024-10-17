function sfm_rawRotationsBuildPairwiseMatrix_test
[E,Rij,RiRef]=sfm_utilityRotationDataset();
G=sfm_rawRotationsBuildPairwiseMatrix(Rij,E);
GRef=sfm_rawRotationsBuildPairwiseMatrixFromR(RiRef,E);
disp(max(abs(vec(G-GRef))))

[~,S,V] = svd(G);
disp(mean(rotationProcrustesAlignError(RiRef,devectorize(V(:,1:3)),'left')))

GLaplacian=sfm_rawRotationsBuildPairwiseMatrix(Rij,E,'method','Laplacian');
[~,S,V] = svd(GLaplacian);
disp(mean(rotationProcrustesAlignError(RiRef,devectorize(V(:,end-2:end)),'left')))


function Ri=devectorize(RiVec)
NRotations=size(RiVec,1)/3;
Ri=reshape(RiVec',[3 3 NRotations]);
Ri=rot_proj(Ri);
