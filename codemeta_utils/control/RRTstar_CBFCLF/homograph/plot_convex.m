function plot_convex(Ax,bx)
s = ["r-","b-","g-","c-"];
x= linspace(0,70,70); 
for i=1:4
    y = (-Ax(i,1)*x+bx(i))/Ax(i,2);
    plot(x,y,s(i))
    hold on
end
axis equal
end