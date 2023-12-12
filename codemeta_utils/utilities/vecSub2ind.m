function idxLinear=vecSub2ind(sz,idx)
szIdx=size(idx);
args=mat2cell(idx,szIdx(1),ones(1,szIdx(2)));
idxLinear=sub2ind(sz,args{:});
