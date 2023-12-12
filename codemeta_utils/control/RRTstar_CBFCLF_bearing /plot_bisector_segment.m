function plot_bisector_segment(tree)
X=0:1:200;
Y=0:1:200;
colors=parula(size(tree,2));
colors = jet(size(tree,2));
idx = randperm(size(tree,2));
for n=1:size(tree,2)
    nn = tree(n).position;
    if ~isempty(nn)
%         plot(nn(1,1),nn(2,1),'k','MarkerSize',5)
        A = tree(n).convex.A;
        b = tree(n).convex.b;
        for i=1:size(X,2)
            for j=1:size(Y,2)
                x = [X(i);Y(j)];
                if all(A*x<=b)
                    plot(x(1),x(2),'Color',colors(idx(n),:),'Marker','.')
                    hold on
                end
            end
        end
    end
end
end
