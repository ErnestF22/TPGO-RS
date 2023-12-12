function POCRelativeDensityBandwidth
phi=@(x) exp(-x.^2/(2*1^2));
%phi=@(x) max(0,1-(x/2).^2);
t=linspace(-4,4);
d=3;
amp=d;
scale=d/4;
P1=phi(t/scale);
P2=amp*phi((t-d)/scale);
plot(t,P1,t,P2,t,P1+P2)
