function POCNormCDF
s=2;
N=1000;
x=s*randn(N,1);
F=@(x) 100*erfc(-x/sqrt(2))/2;
% cumDistPerc(x);
% hold on
% ezplot(F)
% hold off
%Fabs=@(x) F(x)-F(-x);
Fabs=@(x) 1-erfc(x/sqrt(2));
% cumDistPerc(abs(x))
% axis([0 5 0 100])
% hold on
% funPlot(@(x) 100*Fabs(x/s),linspace(0,5))
% hold off

normalQuantile=@(F) erfcinv(1-F)*sqrt(2);
disp(normalQuantile(Fabs(1)))


xSorted=sort(abs(x));
% K=round(0.1*N);
% KPerc=K/N;
% xK=xSorted(K);
% disp([KPerc Fabs(xK/s)])
% disp([s scaleEstimate(xSorted,K,N)])

sK=@(K) scaleEstimate(xSorted,K,N);
funPlot(sK,1:0.5*N)

%vQuantile=@(p) p*(1-p)/(N*normal(normalQuantile(p))^2)
ax=axis();
ax(3:4)=[1 3];
axis(ax)
hold on

end
