function POCModifiedFrobeniuousRotationDistance
I=eye(3);
R=rot_randGeodFun(I);

fRef=@(t) norm(I-R(t),'fro')^2/8;
switch 2
    case 1
        f0=@(t) norm(100*I-R(t),'fro')^2;
    case 2
        f0=@(t) norm(I-R(t),'fro')^2+0.23*norm(R(t)-R(t)','fro');
end

f=@(t) (f0(t)-f0(0))/(f0(pi)-f0(0));
t=linspace(0,pi,70);

funPlot(@(t) fRef(t),t)
hold on
funPlot(@(t) f(t),t,'r')
funPlot(@(t) (1-cos(t))/2,t,'kx')
hold off
