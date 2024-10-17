function v=essential_randTangentNormVector(Q)
N=size(Q,3);
v=zeros(6,3,N);
for iN=1:N
    Q1=Q(:,:,iN);
    v1=essential_tangentProj(Q1,randn(6,3));
    v1=v1/sqrt(essential_metric(Q1,v1,v1));
    v(:,:,iN)=v1;
end
