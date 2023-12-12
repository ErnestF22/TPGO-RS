function plot_circle(x,y,r,s)
th = 0:pi/50:2*pi;
xunit = r * cos(th) + x;
yunit = r * sin(th) + y;
fill(xunit, yunit,s,'facealpha',.2);
hold on
end