function z_neighbor=ReluPlexNeighborsBit(z)
for i=1:size(z,2)
    z_neighbor(i,:) = z;
    z_neighbor(i,i) = double(~(z(i)));
end
z_neighbor = [z;z_neighbor];
end