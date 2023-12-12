function tree = find_controller(tree,l1,l2,l3,l4)
% y=[20 45;40 45;25 5;2 15]';
Au = [1 0;0 1;-1 0;0 -1];
bu = 100*[5;5;5;5];
b = 0*ones(1,size(tree,2));
for i = 1:size(tree,2)
    if tree(i).parent ==1
        b(i) = 0;
    end
end
Wb = ones(size(tree(i).CBF.Ah,1),1);
Wl = 1;
Cb = 1;
Cl = 1;
for i=1:size(tree,2)
    if ~isempty(tree(i).parent)
        y = tree(i).L;
        if y == 1
            y = l1;
        elseif y==2
            y = l2;
        elseif y==3
            y = l3;
        else
            y = l4;
        end
        [K] = optFirstorder(Wb,Wl,...
            Cb,Cl,-tree(i).CLF.z,tree(i).CBF.Ah,Au,...
            tree(i).convex.A,tree(i).convex.b,tree(i).CBF.bh,bu,y,tree(i).CLF.xe,b(i));
         tree(i).K = K;
    end
end
end