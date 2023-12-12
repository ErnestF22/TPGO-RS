% For each node in near_inds, this will check if it is cheaper to
% arrive to them from new_node.
% In such a case, this will re-assign the parent of the nodes in
% near_inds to new_node.
function tree = rewire(tree,x_new,near_indeces,obstacle_list,path_resolution)
for inear=1:size(near_indeces,2)
    x_near = tree(near_indeces(inear)).position;
    cflag = check_collision(x_near,x_new,obstacle_list,path_resolution);
    newCost = tree(end).cost+norm(x_new-tree(near_indeces(inear)).position);
    if cflag &&(newCost<tree(near_indeces(inear)).cost)
        tree(near_indeces(inear)).parent = size(tree,2);
        tree = propagate_cost_to_leaves(tree, near_indeces(inear));
    end
end
end