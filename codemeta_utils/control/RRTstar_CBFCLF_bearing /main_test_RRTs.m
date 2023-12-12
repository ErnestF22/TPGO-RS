% this function gives you the simplified version of RRT*
function [path,fatal_samples,T,tree,rnd]= main_test_RRTs(obstacle_list,start,goal,...
    max_itr,lim_x,lim_y,expand_dis,path_resolution,flag)


tree(1).position = start;
tree(1).cost = 0;
tree(1).parent = [];
tree(1).children = [];

figure(1)
% for i=1:size(obstacle_list,2)
%     plot_circle(obstacle_list(1,i),obstacle_list(2,i)...
%         ,obstacle_list(3,i),'b')
% %     plot_circle(obstacle_list(1,i),obstacle_list(2,i)...
% %         ,10,'r')
% end

if flag
    axis equal
    axis ([0 100 0 100])
else
    axis equal
    axis ([43 475 89 590])
end

figure(2)
% for i=1:size(obstacle_list,2)
%     plot_circle(obstacle_list(1,i),obstacle_list(2,i)...
%         ,obstacle_list(3,i),'b')
% end
if flag
    axis equal
    axis ([0 100 0 100])
else
    axis equal
    axis ([43 475 89 590])
end

[path,fatal_samples,tree,rnd] = expanding_tree(tree,goal,max_itr,lim_x,lim_y ...
    ,expand_dis,path_resolution,obstacle_list);

figure(1)
draw_graph(tree, start,'g-')

for i=1:size(tree,2)
    x1 = tree(i).position;
    if ~isempty(x1)
        plot(x1(1),x1(2),'o','MarkerEdge','k','MarkerFace','k','MarkerSize',5)
    end
end
x1 = tree(1).position;
plot(x1(1),x1(2),'s','MarkerEdge','k','MarkerFace','k','MarkerSize',15)


test = 0;
if test
    figure(2)
    tree = PPR2(tree,obstacle_list,path_resolution);
    draw_graph(tree, start,'g-')
    T = keep_middle_leaf(tree);
    T = RemoveCrossing2(T);
    T = keep_middle_leaf(T);
    
    %
    for i=1:5
        T = PPR2(T,obstacle_list,path_resolution);
        T = keep_middle_leaf(T);
        T = RemoveCrossing2(T);
        T = keep_middle_leaf(T);
        T = RemoveCrossing2(T);
        T = keep_middle_leaf(T);
    end
    
    
    plot_simplified_tree(T,'b',1.5,2);
else
    T = [];
end
end

