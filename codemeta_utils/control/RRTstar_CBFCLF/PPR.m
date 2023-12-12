% post processing rewiring function
function tree = PPR(tree,obstacle_list,path_resolution)

for i=2:size(tree,2)
    if i==330
        a=100;
    end
    if ~isempty(tree(i).position)
        loc_node = tree(i).position;
        T = 1;
        while ~isempty(tree(i).parent) && ~isempty(tree(tree(i).parent).parent) && T 
            idx_parent = tree(i).parent;
            idx_parent_parent = tree(idx_parent).parent;
            loc_parent = tree(idx_parent_parent).position;
            if check_collision(loc_node,loc_parent,obstacle_list,path_resolution)
                tree(i).parent = idx_parent_parent;
                tree(i).cost = tree(idx_parent_parent).cost + ...
                    norm(tree(i).position - tree(idx_parent_parent).position);
            else
                T=0;
            end
        end
    end
end
end
