function tree= fitLandmarks(tree)
for i=1:size(tree,2)
    tree(i).L = 1;
%     if ~isempty(tree(i).parent) && ~isempty(tree(i).position)
%         x = (tree(tree(i).parent).position);
%         %  if (0<=x(1) && x(1)<20) && (0<=x(2)&& x(2)<100)
%         if (0<=x(1) && x(1)<130) && (0<=x(2)&& x(2)<110)
%             tree(i).L = 1;
%             % tree(i).L = l4;
%             % elseif (80<=x(1) && x(1)<=100) &&  (0<=x(2) && x(2)<=20)
%         elseif (130<=x(1) && x(1)<=190) &&  (0<=x(2) && x(2)<=70)
%             tree(i).L = 2;
%             %             tree(i).L = l2;
%             
%             %         elseif (20<=x(1) && x(1)<80) &&  (0<=x(2) && x(2)<20)
%         elseif (130<=x(1) && x(1)<190) &&  (70<=x(2) && x(2)<170)
%             tree(i).L = 3;
%             %             tree(i).L = l3;
%         else
%             tree(i).L = 4;
%             %             tree(i).L = l1;
%         end
%         
%     end
end
end