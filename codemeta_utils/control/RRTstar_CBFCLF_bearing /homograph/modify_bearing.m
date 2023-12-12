function yp = modify_bearing(L,y,idx)
T = landmarksRelativeAngles(L);
Lp = y(:,idx);
M = y(2,:)./y(1,:);
b = zeros(1,size(T,2));
for i=2:size(T,2)
    b(i) = Lp(2)-T(i)*Lp(1);
end
yp = y;
for i=2:size(T,2)
    yp(:,i) = [b(i)/(M(i)-T(i));M(i)*b(i)/(M(i)-T(i))];
end
figure(3)
plot(0,0,'g*')
hold on 
r = 1;
th = 0:pi/50:2*pi;
xunit = r * cos(th) ;
yunit = r * sin(th) ;
h = plot(xunit, yunit);
hold on
for i=1:size(T,2)
    plot(y(1,i),y(2,i),'r*')
    hold on 
    plot(yp(1,i),yp(2,i),'b*')
    hold on 
    plot(L(1,i),L(2,i),'c*')
end
end