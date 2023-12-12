clear all
close all
[tree,~] = testData;
% text(x(1,:)+0.5,x(2,:)+0.5,string(a(:)))
tree = voronoi(tree);
plot_bisector_segment(tree)