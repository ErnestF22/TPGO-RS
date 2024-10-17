function d=multidet(A)
NA=size(A,3);
d=zeros(NA,1);
for iA=1:NA
    d(iA)=det(A(:,:,iA));
end
