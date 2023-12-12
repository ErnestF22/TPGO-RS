function tree= fitLandmarks_exp(tree)
for i=1:size(tree,2)
    if ~isempty(tree(i).parent) && ~isempty(tree(i).position)
        x = (tree(tree(i).parent).position);   
        if (0<=x(1) && x(1)<214) && (0<=x(2)&& x(2)<318)
            tree(i).L = 1;           
        elseif (214<=x(1) && x(1)<=390) &&  (0<=x(2) && x(2)<=220)
            tree(i).L = 2;
        elseif (214<=x(1) && x(1)<390) &&  (220<=x(2) && x(2)<420)
            tree(i).L = 3;
        else
            tree(i).L = 4;
        end
        
    end
end
end