function [vertices,ACTVs]=ReluPlexFindVertexQuery(basic,Ar,d)
k = size(Ar,2)-1;
Q = sort(basic,'ascend'); %tried basic set
B = basic;
T = {Ar}; %tableu set
result = zeros(1,size(Ar,2)-1);
result(1,basic(:))=Ar(1:(size(Ar,1)-1),end);
v = result(1:d)-result(d+1:2*d);
vertices = v;
idx = (2*d+1:size(Ar,2)-1);
actv = idx(~ismember(idx,basic));
ACTVs = actv;
i=1;
while i <= size(B,1) % check if visit all corners
    basic = B(i,:);
    base = sort(B(i,:),'ascend'); % reorder

    A = T{i};
    
    for xs = base % for every variable in basic
        pivotRow = find(basic==xs);
        % find pivot column
        for j=1:k
            if any(base==j) % if already in basic
                continue
            elseif (xs<=(2*d) && j>(2*d)) % if try to pivot basic to slack
                continue
            elseif ((j<=2*d) && (ismember(j+d,base) || ismember(j-d,base)))    
                continue
            else
                % find new base
                base_new = base;
                base_new(base_new == xs)=j;
                base_new = sort(base_new,'ascend');
                act_new = idx(~ismember(idx,base_new));
                
%                 if ismember(base_new,Q,'rows')% if already in Q
                if ismember(act_new,ACTVs,'rows')% if already in Q
                    continue
                else % if not in Q
%                   Q = [Q;base_new]; % put in Q
%                   ACTVs = [ACTVs;act_new];
                    pivotColumn = j;
                    if A(pivotRow,pivotColumn)<= 1e-5 % if the pivot value =0, decided by A
                        continue
                    else %pivot
                        A_new = pivot(A,pivotRow,pivotColumn);
                        basic_new = basic;
                        basic_new(pivotRow)=pivotColumn;
                        [basic_new,~,A_new] = DualSimplex(A_new,basic_new,2);
                        if ~isempty(basic_new) % if it is feasible, then find a new vertex
                            B = [B;basic_new];
                            T{end+1} = A_new;
                            result = zeros(1,size(A,2)-1);
                            result(1,basic_new(:))=A_new(1:(size(A_new,1)-1),end);
                            v = result(1:d)-result(d+1:2*d);
                            actv = idx(~ismember(idx,basic_new));
                            if ~ismember(v,vertices,'rows')
                                vertices = [vertices;v];
                                ACTVs =  [ACTVs;actv];
                            end
                        end
                    end
                end
            end
        end
    end
    i = i+1;
end
end













