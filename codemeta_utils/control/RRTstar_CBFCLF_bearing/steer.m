function nodeSteer = steer(from_node,to_node,extend_length)
if norm(from_node-to_node)>extend_length
    nodeSteer = from_node + (extend_length*(to_node-from_node))/norm(to_node-from_node);
else
    nodeSteer = to_node;
end
end