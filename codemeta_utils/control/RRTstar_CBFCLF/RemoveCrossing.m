%remove crossings
function T = RemoveCrossing(T)
i = 0;
while i<size(T,2)
    i=i+1;
    if ~isempty(T(i).position) %&& checkLeaves(T,i)
        j = 0;
        while j<size(T,2)
            j = j+1;
            if ~isempty(T(j).position) && i~=j
                if  ~isempty(T(i).parent) && ~isempty(T(j).parent)
                    p1 = T(i).position;
                    q1 = T(T(i).parent).position;
                    p2 = T(j).position;
                    q2 = T(T(j).parent).position;
                    flag1 = all(isnan(q1));
                    flag2 = all(isnan(q2));
                    if flag1 || flag2
                        a=1000;
                    end
                    if doIntersect(p1,q1,p2,q2)
                        T = modefied_simplifying(T,i,j);
                    end
                end
            end
        end
    end
end
end