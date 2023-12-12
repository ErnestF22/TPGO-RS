function k_idx  = find_nearest_k(tree,x,obstacle)
dist=[];
id=[];
for i=1:size(tree,2)
    if ~isempty(tree(i).position) && ~isempty(tree(i).parent)
        n = tree(i).position;
        p = tree(tree(i).parent).position;
        if check_collision(p,x,obstacle,1) %&& check_collision(n,x,obstacle,1)
            d = pointToLine(x,n,p);
            d = norm(n-x);
            dist = [dist d];
            id = [id i];
        end
    end
end
[~,idx] = min(dist);
k_idx = id(idx);
end