function plotLines(a,b,c,s)
x = linspace(-100,100,200);
y = -(a/b)*x-(c/b);
plot3(x,y,zeros(size(x,2),1),s,'linewidth',3);
end