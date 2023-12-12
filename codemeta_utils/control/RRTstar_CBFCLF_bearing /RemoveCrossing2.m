%remove crossings
function T = RemoveCrossing2(T,obstacle_list,path_resolution)
i = 0;
while i<size(T,2)
    i=i+1;
    if ~isempty(T(i).position) 
        j = 0;
        while j<size(T,2)
            j = j+1;
            if ~isempty(T(j).position) && i~=j
                if  ~isempty(T(i).parent) && ~isempty(T(j).parent)
                    if (T(i).parent~=j && T(j).parent~=i) && (T(j).parent~=T(i).parent)
                        p1 = T(i).position;
                        q1 = T(T(i).parent).position;
                        p2 = T(j).position;
                        q2 = T(T(j).parent).position;
                        if doIntersect(p1,q1,p2,q2)
                            T = modefied_simplifying(T,i,j);
                        end
                    end
                end
            end
        end
    end
end
end
%
%                             parent1 = T(T(i).parent).position;
%                             parent2 = T(T(j).parent).position;
%                             F1 = check_collision(parent1,T(j).position,obstacle_list,path_resolution);
%                             F2 = check_collision(parent2,T(i).position,obstacle_list,path_resolution);
%                             if F1 && F2
%                                 tmp = T(i).parent;
%                                 parent_i_children = T(T(i).parent).children;
%                                 idx_i = find(parent_i_children==i);
%                                 parent_i_children(idx_i)=j;
%                                 T(T(i).parent).children = parent_i_children;
%                                 
%                                 parent_j_children = T(T(j).parent).children;
%                                 idx_j = find(parent_j_children==j);
%                                 parent_j_children(idx_j)=i;
%                                 T(T(j).parent).children = parent_j_children;
%                                 T(i).parent = T(j).parent;
%                                 T(j).parent = tmp;
%                                 
%                             end

