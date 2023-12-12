function S = find_siblings(parent,idx,node)
S = [];
for i=1:size(node,2)
    if ~(isempty(node(i).parent))
        if i ~= idx
            if (node(i).parent == parent)
                S=[S;i];
            end
        end
    end
end
p = node(parent).parent;
S=[S;p];
end