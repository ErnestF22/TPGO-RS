function rot3_normLogDiff_test

[Rt,dRt,R0,dR0,v]=rot_randGeodFun(rot_randn(eye(3)));

LogRt=@(t) logrot(Rt(t));
theta=@(t) norm(LogRt(t));
u=@(t) LogRt(t)/theta(t);

Du=@(t) rot3_normLogDiff(eye(3),Rt(t));

du=@(t) Du(t)*v;

eigDu=@(t) eig(Du(t));
eigDuTheory=@(t) 0.5*[0;1j+cot(theta(t)/2);-1j+cot(theta(t)/2)];

figure(1)
%check_der(thetat,dthetat,'angle')
check_der(u,du)

figure(2)
subplot(2,1,1)
plotfun(@(t) real(eigDu(t)),'angle','r')
hold on
plotfun(@(t) real(eigDuTheory(t)),'angle','gx')
hold off
title('Real part')
subplot(2,1,2)
plotfun(@(t) imag(eigDu(t)),'angle','r')
hold on
plotfun(@(t) imag(eigDuTheory(t)),'angle','gx')
hold off
title('Imaginary part')
legend('Computed','Expected from theory')

