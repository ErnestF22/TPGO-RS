%Given the tree and samples in obstacles, this function findes the CBF
%constraints for each node of the tree and save them in tree.CBF 
function tree = CBF_constraint(tree,obs)
for i=1:size(tree,2)
    if ~isempty(tree(i).parent)
        n = tree(i).position;
        p = tree(tree(i).parent).position;
        o = obs(:,tree(i).nearest);
        [A,b] = fitParallel(n,p,o);
%         [A,b] = fit_line_with_two_point(n,o);
        [Ax,bx] = make_set_convex(A,b,p);
        tree(i).CBF.Ah = -Ax;
        tree(i).CBF.bh = bx;

        
        [A,b] = fitPerpendicular(n,p);
        [Ax,bx]=make_set_convex(A,b,p);
        tree(i).CBF.Ah = [tree(i).CBF.Ah ;-Ax];
        tree(i).CBF.bh = [tree(i).CBF.bh ;bx];
        
    end
end
end