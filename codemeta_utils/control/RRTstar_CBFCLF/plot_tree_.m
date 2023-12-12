function plot_tree_(tree)
for i=1:size(tree,2)
    if ~isempty(tree(i).position)
        plot(tree(i).position(1),tree(i).position(2),'ro')
        text(tree(i).position(1)+0.5,tree(i).position(2)+0.5,string(i))
        hold on
        if ~(isempty(tree(i).parent))
            plot([tree(i).position(1) tree(tree(i).parent).position(1)],...
                [tree(i).position(2) tree(tree(i).parent).position(2)],...
                'b-','LineWidth',1)
            hold on
        end
    end
end
% plot(obs(1,:),obs(2,:),'bx','MarkerSize',15)
% hold on
% axis equal
% axis([0 60 5 60])
end