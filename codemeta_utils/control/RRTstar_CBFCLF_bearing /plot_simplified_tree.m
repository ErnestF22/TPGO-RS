function t = plot_simplified_tree(T,s,L,n)
t(1).position = T(1).position;
t(1).parent = [];
for i=2:size(T,2)
    if ~isempty(T(i).position)
%         if  ~isempty(T(T(i).parent).position)
            x1 = T(i).position;
            x2 = T(T(i).parent).position;
            
            t(end+1).position =  T(i).position;
            t(end).parent =  T(i).parent;
            figure(n)
            plot([x1(1) x2(1)],[x1(2) x2(2)],s,'linewidth',L);
            %'Color', [0.6 0.6 0.6]
            hold on
%         end
    end
end
end