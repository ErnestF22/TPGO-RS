function locSegGatherClicks
axis([-5 5 -5 5])
grid on
[u,v]=ginput();
x=[u v]';
display(x)
