% CLF constrains:
function tree = CLF_constraint(tree)
for iNode=1:size(tree,2)
    if ~isempty(tree(iNode).parent)
        p = tree(iNode).parent;
        z = tree(iNode).position-tree(p).position;
        tree(iNode).CLF.z = z/norm(z,2);
        tree(iNode).CLF.xe = (tree(p).position);         
    end
end
end