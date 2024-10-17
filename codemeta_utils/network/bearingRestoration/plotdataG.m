clc
clear
close all

load('DATAGreedy1050.mat')

figure
hold on;
%color=['y','m','c','r','g','b','k'];
color=jet(9);
N=10:5:50;
for Nidx=1:length(N)
    n=N(Nidx);
    IntervalI=n-1;
    IntervalF=floor(0.6*(n-1)*(n-2)/2);
    IntervalLength=IntervalF-IntervalI+1;
    
    TaveG=Tave{Nidx};
    
    plot((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'color',color(Nidx,:),'DisplayName',['n=',num2str(n), ' G']);
end

lgd=legend('show','Location','eastoutside');
lgd.FontSize = 8;
xlabel('D');
ylabel('Running Time (s)');