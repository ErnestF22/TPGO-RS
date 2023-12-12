%Computes the cheapest point to new_node contained in the list
%near_inds and set such a node as the parent of new_node.
function [nearest_idx,cost_min] = choose_parent(x_new,nearest_idx,cost_min,near_indeces, obstacle_list,tree,path_resolution)

if ~isempty(near_indeces)
    for inear=1:size(near_indeces,2)
        x_near = tree(near_indeces(inear)).position;
        cflag = check_collision(x_near,x_new,obstacle_list,path_resolution);
        newCost = tree(near_indeces(inear)).cost + ...
            norm(x_new-tree(near_indeces(inear)).position);
        if cflag && (newCost<cost_min)
            cost_min = newCost;
            nearest_idx = near_indeces(inear);
        end
    end
end

end