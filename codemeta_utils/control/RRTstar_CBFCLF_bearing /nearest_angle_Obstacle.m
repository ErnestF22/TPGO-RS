function tree = nearest_angle_Obstacle(tree, obs)
for iNode=1:size(tree,2)
    tree(iNode).nearest_ang = [];
    n = tree(iNode).position;
    if ~isempty(tree(iNode).parent)
        p = tree(tree(iNode).parent).position;
        d=[];
        for j=1:size(obs,2)
            dis = pointToLine(obs(:,j),n,p);
            d = [d,dis];
        end
        if ~isempty(d)
           [ang_l,ang_r] = compute_angle_obs (n,p,obs,d);
           tree(iNode).nearest_ang = [ang_l,ang_r];
        end
    end
end
end