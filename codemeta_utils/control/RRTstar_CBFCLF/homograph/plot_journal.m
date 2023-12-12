function plot_journal(M,K,start,flag_loop,fig)
figure(fig)
%plot env:
for t = 1:size(M,2)
    y = M(t).y;
    for i=1:size(y,2)
        plot(y(1,i),y(2,i),'b.')
        hold on
    end
    plot( y(1,:), y(2,:), '.-b');
    hold on
    plot( [y(1,end) y(1,1)],[y(2,end) y(2,1)], '.-b');
%     plot_controllers(M(t).Ax,M(t).bx,K(t).K,K(t).added,M(t).L,5,0,0,1)
end
S = [1 .6 .6;.5 .5 1];
for t=1:size(M,2)
    L = M(t).L;
    st = S(fix((t-1)/3)+1,:);
    for i=1:size(L,2)
        plot(L(1,i),L(2,i),'s','MarkerSize',10,...
            'MarkerEdgeColor',st,...
            'MarkerFaceColor',st)%st(fix((i-1)/3)+1)
        hold on
    end
end
plot( start(1), start(2), 'o','MarkerSize',6,...
    'MarkerEdgeColor','black',...
    'MarkerFaceColor','black');

navigate_journsl(start,K,M,'m.',flag_loop,fig);

hold on 
% plot( 28, 25, 'd','MarkerSize',6,...
%     'MarkerEdgeColor','black',...
%     'MarkerFaceColor','black');

% set(gcf,'Renderer','Painters');
% setFigFontSize(16)
hold off
end