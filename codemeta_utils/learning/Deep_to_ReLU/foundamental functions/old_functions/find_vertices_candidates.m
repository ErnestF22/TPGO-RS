function Q = find_vertices_candidates(A,basic)
Q = struct([]);
idx = 1:size(A,2)-1;
idxNonBasic = idx(~ismember(idx,basic));
idxRows = 1:size(A,1)-1;
c = 1;
for i=1:size(idxNonBasic,2)
    pivotCol = idxNonBasic(i);
    rowCandidate = idxRows(A(:,pivotCol)>10e-5);
    if ~isempty(rowCandidate)
        for j=1:size(rowCandidate,2)
            Q(c).table = A;
            Q(c).pivotCol = pivotCol;
            Q(c).pivotRow = rowCandidate(j);
            Q(c).basic = basic;
            c=c+1;
        end
    end
end
end