% given the set of vertices, this function plots them

function plot_vertices(vertices,Ab_set,L,s)
for i=1:size(vertices,1)
    [vertices(i,3),~]=ReluPlexNetworkOutput(vertices(i,1:2)',Ab_set,L);
    plotPoints(vertices(i,:)',s,'Markersize',10);
    hold on;
end