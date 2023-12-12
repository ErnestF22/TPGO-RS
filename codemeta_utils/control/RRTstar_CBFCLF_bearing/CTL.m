% cutting the leaves 
function T = CTL(tree)
T = tree;
for i =1: size(tree,2)
    flag = 1;
    for j =1: size(tree,2)
        if tree(j).parent == i
            flag =0;
        end
    end
    if flag
        T(i).position = [];
        T(i).parent = [];
    end
end
end