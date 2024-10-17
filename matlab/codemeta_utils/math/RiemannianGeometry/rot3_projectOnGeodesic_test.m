function rot3_projectOnGeodesic_test
%See geodesicMinimizationSO3.tex for derivations

R=rot_randn();
Rz0=rot_randn();
vRz0=rot_randTangentNormVector(Rz0);
Rz=@(t) rot_exp(Rz0,t*vRz0);
f=@(t) rot_dist(R,Rz(t));

RProjAlgebraic=rot3_projectOnGeodesicAlgebraic(R,Rz0,vRz0);
fProjAlgebraic=rot_dist(R,RProjAlgebraic);

[RProjTrigonometric,tProj]=rot3_projectOnGeodesicTrigonometric(R,Rz0,vRz0);
fProjTrigonometric=f(tProj);

funPlot(f,'angle','k')
axis([-pi pi 0 pi])
hold on
plot([-pi pi],[fProjAlgebraic fProjAlgebraic])
plot([tProj tProj], [0 fProjTrigonometric],'g')
hold off

legend('Distance function','Algebraic solution','Trigonometric solution')

disp('[RProjAlgebraic RProjTrigonometric]')
disp([RProjAlgebraic RProjTrigonometric])

disp('[rot_dist(R, RProjAlgebraic) rot_dist(R, RProjTrigonometric)]')
disp([rot_dist(R, RProjAlgebraic) rot_dist(R, RProjTrigonometric)])
