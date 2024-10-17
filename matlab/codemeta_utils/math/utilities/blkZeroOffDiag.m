%Put to zero all the entries off the block-diagonal of a matrix
%function D=blkZeroOffDiag(A,d)
%Inputs
%   A   Initial matrix
%   d   Size of the blocks (must divide size of A exactly)
function D=blkZeroOffDiag(A,d)
NA=min(size(A));
NBlocks=round(NA/d);
if NBlocks*d~=NA
    error('Size of matrix must be exact multiple of the size of blocks')
end

idxBlock=reshape(1:NA,d,NBlocks);
D=zeros(size(A));
for iBlock=1:NBlocks
    D(idxBlock(:,iBlock),idxBlock(:,iBlock))=A(idxBlock(:,iBlock),idxBlock(:,iBlock));
end