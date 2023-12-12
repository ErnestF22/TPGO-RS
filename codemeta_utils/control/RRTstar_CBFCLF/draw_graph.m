% this function plot the tree with style 's'
function draw_graph(tree, start,s)
for i=1:size(tree,2)
    x1 = tree(i).position;
    if ~isempty(tree(i).parent)
        parent = tree(i).parent;
        x2 = tree(parent).position;
        plot([x1(1) x2(1)],[x1(2) x2(2)],s,'LineWidth',1)
       
        
    end
end
plot(start(1),start(2),'*','MarkerSize',10,...
    'MarkerEdgeColor','black')
% plot(goal(1),goal(2),'*','MarkerSize',10,...
%     'MarkerEdgeColor','blue')
end