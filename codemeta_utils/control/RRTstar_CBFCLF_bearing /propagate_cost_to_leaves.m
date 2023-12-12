function tree = propagate_cost_to_leaves(tree, parent_node)
for i = 1:size(tree,2)
    if tree(i).parent == parent_node
        tree(i).cost = calc_new_cost(tree(parent_node),tree(i).position);
        tree = propagate_cost_to_leaves(tree, i);
    end
end
end