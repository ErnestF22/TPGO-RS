function mDiag=multidiag(m)
nbSlices=size(m,3);
dVector=numel(m(:,:,1));
mDiag=zeros(dVector,dVector,nbSlices);
for iSlice=1:nbSlices
    mDiag(:,:,iSlice)=diag(m(:,:,iSlice));
end
end %file function