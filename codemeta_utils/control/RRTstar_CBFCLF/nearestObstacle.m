function tree = nearestObstacle(tree, obs)
for iNode=1:size(tree,2)
    tree(iNode).nearest = [];
    n = tree(iNode).position;
    if ~isempty(tree(iNode).parent)
        p = tree(tree(iNode).parent).position;
        d=[];
        for j=1:size(obs,2)
            dis = pointToLine(obs(:,j),n,p);
            d = [d,dis];
        end
        if ~isempty(d)
            [~,idx] = min(d);
            tree(iNode).nearest = idx;
            tree(iNode).dis = d(idx);
        end
    end
end
end