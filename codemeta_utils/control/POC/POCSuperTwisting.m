function POCSuperTwisting
xr=-1;
sqrtSign=@(x) sqrt(abs(x)).*sign(x);
s=@(t,x) 1*[-1*abs(x(1,:));-0*ones(size(x(2,:)))];
switch 2
    case 1
        k1=1;
        k3=1;
        dx=@(t,x) [...
            -k1*sqrtSign(x(1,:)-xr)+x(2,:);
            -k3*sign(x(1,:)-xr)
            ]+s(t,x);
    case 2
        k1=6;
        k3=1;
        dx=@(t,x) [...
            -k1*(x(1,:)-xr)+x(2,:);
            -k3*sign(x(1,:)-xr)
            ]+s(t,x);
    case 3
        k1=6;
        k3=1;
        dx=@(t,x) [...
            -k1*sqrtSign(x(1,:));
            zeros(size(x(2,:)))
            ];
        
end

TFinal=5;
x0=[1;0];

optsOde=odeset('OutputFcn',@odeplot,'MaxStep',0.01);
switch 2
    case 1
        [t,x]=ode45(dx,[0 TFinal],x0,optsOde);
    case 2
        [t,x]=odeEuler(dx,[0 TFinal],x0,optsOde);
end
figure(1)
subplot(2,1,1)
plot(t,x)
legend('x_1','x_2')
grid on
subplot(2,1,2)
plot(t,dx(0,x'))
legend('dx_1','dx_2')
grid on
