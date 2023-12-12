%find nearest node from all nodes of the tree
function min_indx=get_nearest_node_index(tree,rnd_node)
for i =1:size(tree,2)
    dis(i) = norm(tree(i).position-rnd_node);
end
[~,min_indx]=min(dis);
end