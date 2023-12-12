function path = generate_path(tree, goal_ind,goal)
path = [goal(1);goal(2)];
path = [path tree(goal_ind).position];
node = tree(goal_ind);
while ~isempty(node.parent)
    path = [path node.position];
    node = tree(node.parent);
end
path = [path node.position];
end