% this function compute the cost of node to_node:
function d = calc_new_cost(from_node, to_node)
[d,~] = calc_distance_and_angle(from_node.position, to_node);
d = d+from_node.cost;
end