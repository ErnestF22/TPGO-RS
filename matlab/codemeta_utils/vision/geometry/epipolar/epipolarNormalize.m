%Divide E by the mean of its first two singular values
%function E=epipolarNormalize(E)
function E=epipolarNormalize(E)
NE=size(E,3);
for iE=1:NE
    s=svd(E(:,:,iE));
    E(:,:,iE)=E(:,:,iE)/mean(s(1:2));
end
