function [V,R] = find_vertices(A,basic)
nVariables = size(A,2)-1;
V = result(A,basic);
variables = 1:1:nVariables;
nonbasic = variables(~ismember(variables,basic));
nonbasicPivot = nonbasic(A(end,nonbasic)==0);
Q(1).rc = query(A,nonbasicPivot);
Q(1).tableau = A;
Q(1).basic = basic;
c = 1;
h=0;
while h<=size(Q,2)
    P = Q(1).rc;
    for k=1:size(P,1)
        [A,basic] = pivoting(Q(1).tableau,P(k,:),Q(1).basic);
        if ~isempty(basic)
            c = c+1;
            new_v = result(A,basic);
            if ~any(ismember(V,new_v,'rows'))
                V = [V;new_v];
                R(c).T = A;
                R(c).V = new_v;
                variables = 1:1:nVariables;
                nonbasic = variables(~ismember(variables,basic));
                nonbasicPivot = nonbasic(A(end,nonbasic)==0);
                Q(c).rc = query(A,nonbasicPivot);
                Q(c).tableau = A;
                Q(c).basic = basic;
            end
        end
    end
    Q(1) = [];
    h=h+1;
end
end