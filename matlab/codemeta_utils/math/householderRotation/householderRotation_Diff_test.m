function householderRotation_Diff_test
mode=2;
[x1,~,~,dx1]=real_randGeodFun(randn(3,1));
switch mode
    case 1
        x2=@(t) 3;
        dx2=[];
    case 2
        [x2,~,~,dx2]=real_randGeodFun(randn(3,1));        
end
check_der(@(t) funDer(x1(t),x2(t),dx1,dx2), 'function')
        

function [R0,dR0]=funDer(x1,x2,dx1,dx2)
R0=householderRotation(x1,x2);
dR0Vec=householderRotation_Diff(x1,x2,dx1,dx2);
dR0=R0*hat(dR0Vec);
