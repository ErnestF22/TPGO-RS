function subM=subMatrix(M,idx1,idx2)
if ~exist('idx2','var')
    idx2=1;
end
subM=M(idx1,idx2);
