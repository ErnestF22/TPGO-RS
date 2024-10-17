function cumDistBoxPerc(data)
cumDistPerc(data)
hAxes1=gca;
colorOrder=get(hAxes1,'ColorOrder');
colorOrder=colorOrder(1:size(data,2),:);

%overimpose box plots
pos=get(gca,'Position');
hAxes2=axes();
boxplot(gca,fliplr(data),'orientation','horizontal','colors',flipud(colorOrder))
pos(2)=1.5*pos(2);
pos(4)=0.8*pos(4);
set(hAxes2,'Position',pos)
set(hAxes2,'YTickLabel',{' '})
%set(hAxes1,'Visible','off')
set(hAxes2,'Visible','off')
set(hAxes1,'YTick',0:20:100,'YGrid','on','YMinorGrid','on','YMinorTick','on','XMinorTick','on','YLim',[0 100])
set(gcf,'CurrentAxes',hAxes1);
