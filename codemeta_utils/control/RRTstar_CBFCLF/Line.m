function path = Line(from_node,to_node,r)
[d,theta] = calc_distance_and_angle(from_node, to_node);
new_node = from_node;
n_expand = floor(d/r);
path = new_node;
for i=1:n_expand
    new_node(1) = new_node(1)+r*cosd(theta);
    new_node(2) = new_node(2)+r*sind(theta);
    path = [path new_node];
end
[d,~] = calc_distance_and_angle(new_node, to_node);
if d<=r
    path = [path to_node];
end
end