function POCNewColormap

N=1024;
I=ones(round(N/5),1)*(1:N);
%c=jet(N);
c=rbg(N);

subplot(2,1,1)
imagesc(I)
colormap(c)
subplot(2,1,2)
plot(c)
hold on
plot([N/2 N/2],[0 1],'k')
hold off
disp(c([1 end],:))
