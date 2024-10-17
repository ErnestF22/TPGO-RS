function G=G2GNormalized(G)
N=size(G,3);
for iN=1:N
    G(1:3,4,iN)=cnormalize(G(1:3,4,iN));
end
