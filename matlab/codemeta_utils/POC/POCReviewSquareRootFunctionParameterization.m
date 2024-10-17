function POCReviewSquareRootFunctionParameterization
flagPlotGraph=false;
%%
syms t

c1=[cos(2*pi*t);sin(2*pi*t)];
c1b=[2*cos(2*pi*t);sin(2*pi*t)];

tau=(((t+1).^2)-1)/2;

c2=subs(c1,t,tau);
c2b=subs(c2,t,tau);

if flagPlotGraph
    ts=linspace(0,1,100);
    plotPoints(double(subs(c1,t,ts)),'.')
    hold on
    plotPoints(double(subs(c2,t,ts)),'x')
    hold off
end

%%
srvf=@(c) diff(c)/(diff(sqrt(norm(c)))+eps);
b1=srvf(c1);
b2=srvf(c2);

ts=linspace(0,1,100);
plotPoints(double(subs(b1,t,ts)),'.')
hold on
plotPoints(double(subs(b2,t,ts)),'x')
hold off


%%