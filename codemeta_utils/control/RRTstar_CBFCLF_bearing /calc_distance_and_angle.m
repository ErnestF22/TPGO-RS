%Given two points from_node and to_node, this function compute the distance
%between two nodes and finds the angle between the x-axis and vector from
%from_node to to_nodes
function [d, theta] = calc_distance_and_angle(from_node, to_node)
v = to_node - from_node;
d = norm(v,2);
%theta = atan(v(2)/v(1));
theta = atan2d(v(2),v(1)) + 360*(v(2)<0);
end
