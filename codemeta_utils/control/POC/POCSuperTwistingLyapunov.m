function POCSuperTwistingLyapunov
sqrtSign=@(x) sqrt(abs(x)).*sign(x);

k1=1;
k3=1;
V=@(x) 2*k3*abs(x(1))+x(2)^2/2 ...
    +(k1*sqrtSign(x(1))-x(2))^2;

x=linspace(-1,1);
[xx1,xx2]=meshgrid(x,x);
xx=cat(3,xx1,xx2);
Vxx=evalfunVec(V,xx);
dx=@(t,x) [...
    -k1*sqrtSign(x(1,:))+x(2,:);
    -k3*sign(x(1,:))
    ];

switch 1
    case 1
        contour(xx1,xx2,Vxx,50)
        hold on
        plotfield(@(x) dx(0,x),x(1:7:end))
        hold off
    case 2
        surf(xx1,xx2,Vxx,'EdgeColor','none')
        zlabel('V')
end
xlabel('x1')
ylabel('x2')
