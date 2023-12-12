%define a ball centered on new_node
%return all nodes of the tree that are inside the ball
function near_inds = find_near_nodes(tree,new_node,r)
near_inds=[];
for i=1:size(tree,2)
    d = norm(tree(i).position-new_node);
    if d<=r
        near_inds = [near_inds i];
%         plot(tree(i).position(1),tree(i).position(2),'m*')
%         hold on
    end
end
end