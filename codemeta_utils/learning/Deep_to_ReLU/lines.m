function lines(A1,A2,A3,b1,b2,b3)

plotLines(A1(1,1),A1(1,2),b1(1),'r')
hold on
plotLines(A1(2,1),A1(2,2),b1(2),'b')
hold on

a = A2*A1;
b = A2*b1+b2;
plotLines(a(1,1),a(1,2),b(1),'g')
hold on
plotLines(a(2,1),a(2,2),b(2),'y')
hold on


a = A3*A2*A1;
b = A3*A2*b1+A3*b2+b3;
plotLines(a(1,1),a(1,2),b(1),'c')

end
