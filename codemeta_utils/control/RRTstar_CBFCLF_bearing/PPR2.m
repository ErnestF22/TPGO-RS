% post processing rewiring function
function tree = PPR2(tree,obstacle_list,path_resolution)

Q = tree(1).children;
visited = [];
while ~isempty(Q)
    
    idx = Q(1);
    q = tree(idx);
    Q(1) = [];
    visited = [visited,idx];
    if ~isempty(q.parent) && ~isempty(q.children)
        Child = q.children;
        parent_idx = q.parent;
        parent_pos = tree(parent_idx).position;
        
        for i=1:size(Child,2)
            child_pos = tree(Child(i)).position;
            if ~any(visited(:) == Child(i))
                if ~isempty(child_pos)
                    if check_collision(child_pos,parent_pos,obstacle_list,path_resolution)
                        tmp = tree(Child(i)).parent;
                        
                        parent_i_children = tree(tmp).children;
                        idx_i = find(parent_i_children==Child(i));
                        parent_i_children(idx_i)=[];
                        tree(tmp).children = parent_i_children;
                        
                        tree(Child(i)).parent = parent_idx;
                        tree(parent_idx).children = [tree(parent_idx).children,Child(i)];
                        tree(Child(i)).cost = tree(parent_idx).cost + ...
                            norm(child_pos - parent_pos);
                    end
                end
                Q = [Q,Child(i)];
            end
        end
    end
end
end
