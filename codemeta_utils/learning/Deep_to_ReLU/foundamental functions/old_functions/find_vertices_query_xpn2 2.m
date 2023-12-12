function [ACT, vertices]=find_vertices_query_xpn2(basic,Ar,d)
% the only thing changed is dual_simplex function(row 41), all the other code
% are exactly the same
qCandidates(1).q = find_vertices_candidates(Ar,basic);
result = zeros(1,size(Ar,2)-1);
result(1,basic(:))=Ar(1:(size(Ar,1)-1),end);
vertices = result(1:d/2)-result(d/2+1:d);
idx = (d+1:size(Ar,2)-1);
a = idx(~ismember(idx,basic));
ACT = sort(a);

while ~isempty(qCandidates) % check if visit all corners
    Q = qCandidates(1).q;
    qCandidates(1) = [];
    for i=1:size(Q,2)
        basic_new = Q(i).basic;
        A = Q(i).table;
        pivotRow = Q(i).pivotRow;
        pivotColumn = Q(i).pivotCol;
        P = eye(size(A,1));
        [A_new,~] = pivotAP(A,P,pivotRow,pivotColumn);
        basic_new(pivotRow) = pivotColumn;
        [basic_new,~,A_new] = dual_simplex_xpn(A_new,basic_new,d);
        a_new = idx(~ismember(idx,basic_new));
        if size(a_new,2)== d/2
%             if ~ismember(sort(a_new),ACT,'rows')
            if sum(sum(ACT-sort(a_new),2)~=0)
                result = zeros(1,size(A_new,2)-1);
                result(1,basic_new(:))=A_new(1:end-1,end);
                v = result(1:d/2)-result(d/2+1:d);
                vertices = [vertices;v];
                ACT = [ACT;sort(a_new)];
                qCandidates(end+1).q = find_vertices_candidates(A_new,basic_new);
            end
        end
    end
end
end