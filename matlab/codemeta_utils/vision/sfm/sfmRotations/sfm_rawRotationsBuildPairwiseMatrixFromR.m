function G=sfm_rawRotationsBuildPairwiseMatrixFromR(R,E)
flagGivenEdges=exist('E','var') && ~isempty(E);
RVec=matStack(permute(R,[2 1 3]));
G=RVec*RVec';
if flagGivenEdges
    NRotations=size(R,3);
    ANeg=~edges2adjmatrix(E,NRotations,'undirected')-eye(NRotations);
    G(logical(kron(ANeg,ones(3))))=0;
end