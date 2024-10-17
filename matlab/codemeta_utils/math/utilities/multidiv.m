%Division of multiple matrices
%function N=multidiv(N,d)
%Loop and compute N(:,:,iN)/d(:,:,iN)
function N=multidiv(N,d)
NN=size(N,3);
for iN=1:NN
    N(:,:,iN)=N(:,:,iN)/d(:,:,iN);
end
