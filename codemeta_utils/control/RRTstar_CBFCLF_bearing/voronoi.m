function tree = voronoi(tree)
for i=1:size(tree,2)
    tree(i).convex.A =[];
    tree(i).convex.b = [];
    AA=[];
    bb=[];
    x = tree(i).position;
    for j=1:size(tree,2)
        y = tree(j).position;
        if (i~=j)
            if ~isempty(x) && ~isempty(y)
                [A,B] = bisector_segment(x,y);
                [Ax,bx] = make_set_convex(A,B,x);
                AA = [AA; Ax];
                bb = [bb; bx];
            end
        end
    end
    if ~isempty(tree(i).parent)
        p = tree(tree(i).parent).position;
        [A,b] = remove_intersection_with_edgeij(AA,bb,x,p);
        tree(i).convex.A = A;
        tree(i).convex.b = b;
        [A,B] = fitPerpendicular(p,x);
        [Ax,bx] = make_set_convex(A,B,x);
        tree(i).convex.A = [tree(i).convex.A; Ax];
        tree(i).convex.b = [tree(i).convex.b; bx];
    end
end
end