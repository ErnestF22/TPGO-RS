function v=lie_randTangentNormVector(lf,y,N)
if ~exist('N','var')
    N=1;
end
v=lf.tangentProj(y,randn([size(y) N]));
for iv=1:size(v,3)
    v(:,:,iv)=v(:,:,iv)/sqrt(lf.metric(y,v(:,:,iv),v(:,:,iv)));
end
