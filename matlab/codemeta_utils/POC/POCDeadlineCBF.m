function POCDeadlineCBF

I=[0,1];
vmax=1;
%ezplot(@(x) dInterval(x,I))
h=@(t,x) t-dInterval(x,I)/vmax;
%fplot(@(x) h(2,x))

B=@(t,x) 1/h(t,x);
fplot(@(x) B(2,x),[-2 3])


%distance from an interval
function d=dInterval(x,I)
d=max(0,x-max(I))+max(0,min(I)-x);

