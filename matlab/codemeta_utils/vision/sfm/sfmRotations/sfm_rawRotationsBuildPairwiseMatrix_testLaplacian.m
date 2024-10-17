function sfm_rawRotationsBuildPairwiseMatrix_testLaplacian

E=[1 2; 2 3];
Rij=rot_randn([],[],2);

Ri=rot_randGeodFun(repmat(eye(3),1,1,3));
RiVec=@(t) reshape(multitransp(Ri(t)),3,[])';

G=sfm_rawRotationsBuildPairwiseMatrix(Rij,E,'method','laplacian');

c=@(t) cost(E,Rij,Ri(t));
cG=@(t) trace(RiVec(t)'*G*RiVec(t));
funCompare(c,cG)

function c=cost(E,Rij,Ri)
NEdges=size(E,1);
c=0;
for iEdge=1:NEdges
    iNode=E(iEdge,1);
    jNode=E(iEdge,2);
    c=c+norm(Ri(:,:,iNode)-Rij(:,:,iEdge)*Ri(:,:,jNode),'fro')^2;
end
