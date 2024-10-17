clc
clear
close all

load('DATASCC.mat')

figure
% hold on;
color=['r','g','b','k'];
for idx=1:3
    n=idx+5;
    IntervalI=n-1;
    IntervalF=floor(0.8*(n-1)*(n-2)/2);
    
    TaveG=Tave1{idx};
    TaveC1=Tave2{idx};
    TaveC2=Tave3{idx};
    LaveG=Lave1{idx};
    LaveC1=Lave2{idx};
    LaveC2=Lave3{idx};
    
    subplot(2,1,1)
    
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveG,'DisplayName',['n=',num2str(n), ' Greedy']);
    hold on
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveC1,'--','color',color(idx),'DisplayName',['n=',num2str(n), ' Comb, First']);
    semilogy((IntervalI:IntervalF)*2/(n*(n-1)),TaveC2,'-*','color',color(idx),'DisplayName',['n=',num2str(n), ' Comb, Best']);

    subplot(2,1,2)
    hold on
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveG,'DisplayName',['n=',num2str(n), ' Greedy']);
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC1,'--','color',color(idx),'DisplayName',['n=',num2str(n), ' Comb, First']);
    plot((IntervalI:IntervalF)*2/(n*(n-1)),LaveC2,'-*','color',color(idx),'DisplayName',['n=',num2str(n), ' Comb, Best']);
end



subplot(2,1,1)
legend('show')
xlabel('D');
ylabel('Running Time (s)');
subplot(2,1,2)
legend('show')
xlabel('D');
ylabel('\sigma_2(M)');