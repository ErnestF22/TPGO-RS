%plot the CBF, the convex set and CLF
close all
x = [2 4;6 2;11 8;2 9]';
j = [8 5];

for i=1:4
    plot([x(1,i) j(1)],[x(2,i) j(2)],'b-')
    hold on 
    plot(x(1,i),x(2,i),'bo')
end
plot(j(1),j(2),'bo')
axis([0 11.5 0 11.5])
axis equal