clc
clear
close all

load('DATASCC789.mat')

figure
% hold on;
color=['r','g','b','k'];
%color=jet(9);
for idx=1:4
    if idx==4
        clear
        idx=4;
        load('DATASCC10.mat')
        Idx=1;
        color=['r','g','b','k'];
    else
        Idx=idx;
    end
    
    n=idx+6;
    IntervalI=n-1;
    IntervalF=floor(0.8*(n-1)*(n-2)/2);
    
    TaveG=Tave1{Idx};
    TaveC1=Tave2{Idx};
    TaveC2=Tave3{Idx};
    LaveG=Lave1{Idx};
    LaveC1=Lave2{Idx};
    LaveC2=Lave3{Idx};
    
    subplot(2,1,1)
    
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'color',color(idx),'DisplayName',['n=',num2str(n), ' G']);
    hold on
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveC1,'--','color',color(idx),'DisplayName',['n=',num2str(n), ' C1']);
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveC2,'-*','color',color(idx),'DisplayName',['n=',num2str(n), ' CB']);

    subplot(2,1,2)
    hold on
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveG,'color',color(idx),'DisplayName',['n=',num2str(n), ' G']);
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC1,'--','color',color(idx),'DisplayName',['n=',num2str(n), ' C1']);
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC2,'-*','color',color(idx),'DisplayName',['n=',num2str(n), ' CB']);
end



subplot(2,1,1)
%legend('show','Location','bestoutside')
lgd=legend('show','Location','eastoutside');
lgd.FontSize = 8;
xlabel('D');
ylabel('Running Time (s)');
subplot(2,1,2)
%legend('show','Location','bestoutside')
lgd=legend('show','Location','eastoutside');
lgd.FontSize = 8;
xlabel('D');
ylabel('\sigma_2(M)');